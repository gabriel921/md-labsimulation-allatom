#!/bin/bash -e
# =============================================================================
# run_setup_analysis.sh
#
# Dispatches analysis jobs (index building + one or more analysis scripts)
# across a user-specified pool of machines via SSH + nohup, since analysis
# cannot run on the same machines/partitions as production here — no SLURM
# for this part of the pipeline.
#
# By default only ONE analysis runs per machine at a time — raise
# MAX_JOBS_PER_MACHINE explicitly if you want to double up.
#
# Run this from the system root (same level as 0_prep/, 4_prd/).
# =============================================================================

set -e

# ── USER CONFIGURATION ────────────────────────────────────────────────────────

# ── EXECUTION MODE ────────────────────────────────────────────────────────────
# "ssh"   — original (default) behaviour: dispatch across MACHINES[] via SSH.
#                      Requires passwordless SSH to every machine in the pool.
# "local" — no SSH, no machine pool. Everything runs on this machine instead,
#           LOCAL_MAX_JOBS job(s) at a time (set LOCAL_MAX_JOBS=1 to run
#           strictly one analysis after another). Use this if you don't have
#           access to other machines. MACHINES[]/SSH_USER are ignored.
EXECUTION_MODE=ssh
LOCAL_MAX_JOBS=1                                      # concurrency when EXECUTION_MODE=local (1 = sequential)

# SSH_USER/MACHINES/MAX_JOBS_PER_MACHINE only apply when EXECUTION_MODE=ssh.
#
# NOTE ON MACHINE NAMES: "bioXXX" below are hostnames specific to this
# author's lab/cluster (bio209, bio210, ...) — they will not exist on your
# network. Replace MACHINES[] with your own reachable hostnames, or just use
# EXECUTION_MODE=local if you don't have a pool of machines to dispatch to.
SSH_USER=""                                           # leave empty to use your current user
MACHINES=(bio209 bio210 bio211 bio212 bio213 bio214)  # pool of machines you're allowed to use THIS run
MAX_JOBS_PER_MACHINE=2                                # raise only if you explicitly want to double up

GROM_BIN_DEFAULT="/gromacs/gromacs-2023.5/bin/gmx"    # gmx binary path used for indexing/analysis

# Per-analysis override, keyed by ANALYSES[] subfolder name. Only needed for
# analyses that are pinned to a specific GROMACS version (e.g. do_dssp's
# argument syntax changed after 2021.x) — anything not listed here falls
# back to GROM_BIN_DEFAULT.
declare -A GROM_BIN_OVERRIDE=(
    ["03_DSSP"]="/gromacs/gromacs-2021.2/bin/gmx"
)

REPLICATES=()                                         # e.g. (rep1 rep2 rep3) — leave empty to auto-discover all repN under 4_prd/
FORCE_REINDEX=0                                       # set to 1 to rebuild index_analysis.ndx even if it already exists
DRY_RUN=0                                             # set to 1 to preview machine assignment without dispatching anything
POLL_INTERVAL=60                                      # seconds between checks for finished jobs, once every slot is full

# Full catalog of known PER-REPLICA analyses. Format per entry:
#   "subfolder:script:extra_file1,extra_file2,..."
# (extra files, e.g. plot_rg.py, get copied alongside the main script;
# leave the trailing ":..." off if there are none)
#
# Master copies of these scripts live at $SCRIPTS_DIR/analysis/<subfolder>/
ANALYSES=(
    "01_Rg:script_Rg_replica.sh:plot_rg.py"
    "02_SASA:script_SASA_replica.sh:plot_sasa.py"
    "03_DSSP:script_dssp_replica.sh:plot_dssp.py"
)

# Which of the entries above to actually run THIS invocation.
# Comma-separated list of subfolder names, e.g. "04_Rg,05_SASA".
# Leave unset/empty to run everything in ANALYSES[].
RUN_ANALYSES="${RUN_ANALYSES:-}"

# Catalog of POOLED analyses — these consume ALL replicas as input and run
# ONCE per invocation (not once per replica). Same "subfolder:script:extras"
# format as ANALYSES[] above, but output lands in $PRD_DIR/Analysis/<subfolder>/
# instead of repN/analysis/<subfolder>/, and only one job is ever enqueued
# per entry regardless of how many replicates exist.
#
# Master copies live at $SCRIPTS_DIR/analysis/<subfolder>/, same convention
# as ANALYSES[].
POOLED_ANALYSES=(
    "04_HDBSCAN:script_hdbscan_pooled.sh:hdbscan_pooled.py"
    "05_Interactions:script_intracontacts_pooled.sh:intra_contacts_pooled.py"
)

RUN_POOLED_ANALYSES="${RUN_POOLED_ANALYSES:-}"        # same filtering convention as RUN_ANALYSES

# nowater.tpr/.gro prep — required once by any pooled analysis. Pinned to
# 2021.2 by default (same rationale as DSSP: version-sensitive gmx tooling).
NOWATER_GROM_BIN="/gromacs/gromacs-2021.2/bin/gmx"
FORCE_NOWATER_REGEN=0                                 # set to 1 to force-regenerate nowater.tpr/.gro even if present

SYSTEM_ROOT="$(pwd)"
SCRIPTS_DIR="$SYSTEM_ROOT/input"                      # system root — analysis/ and 00_index/ live under here
PRD_DIR="$SYSTEM_ROOT/4_prd"
POOLED_DIR="$PRD_DIR/Analysis"
# ── END USER CONFIGURATION ───────────────────────────────────────────────────

log() { echo -e "[$(date '+%H:%M:%S')] $*"; }
die() { echo -e "[ERROR] $*" >&2; exit 1; }

case "$EXECUTION_MODE" in
    ssh|local) ;;
    *) die "Unknown EXECUTION_MODE='$EXECUTION_MODE' — must be 'ssh' or 'local'." ;;
esac

# In local mode there is no machine pool: everything runs on this machine,
# LOCAL_MAX_JOBS at a time. Reusing the existing MACHINES[]/MAX_JOBS_PER_MACHINE
# machinery with a single pseudo-machine means the queue/dispatch/poll logic
# below (unchanged) works identically for both modes.
if [ "$EXECUTION_MODE" == "local" ]; then
    MACHINES=(LOCAL)
    MAX_JOBS_PER_MACHINE=$LOCAL_MAX_JOBS
fi

[ -d "$PRD_DIR" ] || die "production directory not found: $PRD_DIR (run this from the system root)"
[ "${#MACHINES[@]}" -gt 0 ] || die "MACHINES[] is empty — specify at least one machine."

ssh_host() {
    if [ -n "$SSH_USER" ]; then echo "${SSH_USER}@$1"; else echo "$1"; fi
}

# ── REMOTE/LOCAL EXECUTION HELPERS ────────────────────────────────────────────
# Every actual point of execution (nowater prep, indexing, dispatch, aliveness
# poll) goes through these three, so EXECUTION_MODE=local runs the exact same
# commands as EXECUTION_MODE=ssh — just via `bash -c` instead of `ssh host`.

# Blocking foreground command. Returns the command's exit status.
remote_exec() {
    local host="$1" cmd="$2"
    if [ "$EXECUTION_MODE" == "local" ]; then
        bash -c "$cmd"
    else
        ssh -o BatchMode=yes "$host" "$cmd"
    fi
}

# Backgrounds job_cmd and prints its PID to stdout. job_cmd should NOT
# include its own trailing "& echo $!" — that's added here per mode.
#
# Local mode wraps job_cmd in its own subshell + redirection before
# backgrounding it. Without this, "job_cmd & echo $!" run via `bash -c`
# would block the PID capture until job_cmd finishes: a `cmd1 && cmd2 &`
# list keeps this bash -c process's own stdout/stderr (the pipe our PID
# capture reads from) open for as long as the backgrounded job runs, even
# though cmd2's own output is separately redirected to the logfile. SSH mode
# is unaffected by this and sends the exact same command string as before.
remote_bg() {
    local host="$1" job_cmd="$2"
    if [ "$EXECUTION_MODE" == "local" ]; then
        bash -c "( $job_cmd ) >/dev/null 2>&1 & echo \$!"
    else
        ssh -o BatchMode=yes "$host" "$job_cmd & echo \$!"
    fi
}

# True (exit 0) if the given PID is still alive on the given machine.
remote_alive() {
    local host="$1" pid="$2"
    if [ "$EXECUTION_MODE" == "local" ]; then
        kill -0 "$pid" 2>/dev/null
    else
        ssh -o BatchMode=yes "$host" "kill -0 $pid" 2>/dev/null
    fi
}

# ── FILTER ANALYSES BY RUN_ANALYSES ───────────────────────────────────────────
if [ -n "$RUN_ANALYSES" ]; then
    IFS=',' read -ra WANTED <<< "$RUN_ANALYSES"
    FILTERED=()
    for entry in "${ANALYSES[@]}"; do
        subfolder="${entry%%:*}"
        for w in "${WANTED[@]}"; do
            if [ "$subfolder" == "$w" ]; then
                FILTERED+=("$entry")
                break
            fi
        done
    done
    [ "${#FILTERED[@]}" -gt 0 ] || die "RUN_ANALYSES=\"$RUN_ANALYSES\" matched nothing in ANALYSES[] (known: $(for a in "${ANALYSES[@]}"; do echo -n "${a%%:*} "; done))"
    ANALYSES=("${FILTERED[@]}")
fi

# ── FILTER POOLED ANALYSES BY RUN_POOLED_ANALYSES ─────────────────────────────
if [ "${#POOLED_ANALYSES[@]}" -gt 0 ] && [ -n "$RUN_POOLED_ANALYSES" ]; then
    IFS=',' read -ra WANTED_POOLED <<< "$RUN_POOLED_ANALYSES"
    FILTERED_POOLED=()
    for entry in "${POOLED_ANALYSES[@]}"; do
        subfolder="${entry%%:*}"
        for w in "${WANTED_POOLED[@]}"; do
            if [ "$subfolder" == "$w" ]; then
                FILTERED_POOLED+=("$entry")
                break
            fi
        done
    done
    [ "${#FILTERED_POOLED[@]}" -gt 0 ] || die "RUN_POOLED_ANALYSES=\"$RUN_POOLED_ANALYSES\" matched nothing in POOLED_ANALYSES[]"
    POOLED_ANALYSES=("${FILTERED_POOLED[@]}")
fi

# ── DISCOVER REPLICATES ───────────────────────────────────────────────────────
if [ "${#REPLICATES[@]}" -eq 0 ]; then
    mapfile -t REPLICATES < <(find "$PRD_DIR" -maxdepth 1 -mindepth 1 -type d -name 'rep*' -exec basename {} \; | sort)
fi
[ "${#REPLICATES[@]}" -gt 0 ] || die "no replicate folders found under $PRD_DIR"
log "Replicates        : ${REPLICATES[*]}"
log "Analyses this run : $(for a in "${ANALYSES[@]}"; do echo -n "${a%%:*} "; done)"
if [ "${#POOLED_ANALYSES[@]}" -gt 0 ]; then
    log "Pooled analyses   : $(for a in "${POOLED_ANALYSES[@]}"; do echo -n "${a%%:*} "; done)"
fi
echo

# ── NOWATER PREP (blocking, once per run, only if a pooled analysis needs it) ─
if [ "${#POOLED_ANALYSES[@]}" -gt 0 ]; then
    NOWATER_MACHINE="${MACHINES[0]}"
    nowater_host=$(ssh_host "$NOWATER_MACHINE")

    log "=== Pooled analyses requested — preparing nowater.tpr/.gro (blocking, on $NOWATER_MACHINE) ==="
    mkdir -p "$POOLED_DIR"
    cp -f "$SCRIPTS_DIR/analysis/prep_nowater.sh" "$POOLED_DIR/"

    nowater_cmd="cd '$POOLED_DIR' && PRD_DIR='$PRD_DIR' GROM='$NOWATER_GROM_BIN' FORCE_REGEN='$FORCE_NOWATER_REGEN' bash prep_nowater.sh"
    if [ "$DRY_RUN" -eq 1 ]; then
        log "  [DRY RUN] would run ($EXECUTION_MODE, $NOWATER_MACHINE): prep_nowater.sh (PRD_DIR='$PRD_DIR', GROM='$NOWATER_GROM_BIN')"
    else
        if ! remote_exec "$nowater_host" "$nowater_cmd"; then
            die "nowater prep failed on $NOWATER_MACHINE — aborting (pooled analyses depend on it)."
        fi
        log "  nowater prep OK."
    fi
    echo
fi

# ── BUILD JOB LIST (replica x analysis, PLUS one job per pooled analysis) ─────
# JOB_TYPE tracks which branch dispatch_job() takes for each queued job:
#   "replica" -> per-replica job (repN/analysis/<subfolder>/), indexed first
#   "pooled"  -> pooled job (4_prd/Analysis/<subfolder>/), no per-job indexing
#                (relies on the nowater prep + rep1's index already existing)
JOB_REPLICAS=()
JOB_ANALYSES=()
JOB_TYPE=()
for rep in "${REPLICATES[@]}"; do
    for entry in "${ANALYSES[@]}"; do
        JOB_REPLICAS+=("$rep")
        JOB_ANALYSES+=("$entry")
        JOB_TYPE+=("replica")
    done
done
for entry in "${POOLED_ANALYSES[@]}"; do
    JOB_REPLICAS+=("POOLED")
    JOB_ANALYSES+=("$entry")
    JOB_TYPE+=("pooled")
done

NUM_JOBS=${#JOB_REPLICAS[@]}
NUM_MACHINES=${#MACHINES[@]}
CAPACITY=$(( NUM_MACHINES * MAX_JOBS_PER_MACHINE ))

log "Jobs to dispatch  : $NUM_JOBS  ($(( ${#REPLICATES[@]} * ${#ANALYSES[@]} )) per-replica + ${#POOLED_ANALYSES[@]} pooled)"
log "Machine capacity  : $CAPACITY concurrent  ($NUM_MACHINES machines x $MAX_JOBS_PER_MACHINE per machine)"
if [ "$NUM_JOBS" -gt "$CAPACITY" ]; then
    log "  Not enough capacity to run all jobs at once — the remainder will queue and backfill as slots free up (poll every ${POLL_INTERVAL}s)."
fi
echo

# ── MANIFEST ───────────────────────────────────────────────────────────────────
MANIFEST="analysis_dispatch_$(date +%Y%m%d_%H%M%S).log"
echo "# replica|analysis|machine|remote_dir|pid|logfile" > "$MANIFEST"

# ── DISPATCH (queue-aware: fills every open slot, polls, backfills) ────────────
QUEUE_IDX=0

declare -A MACHINE_LOAD
for m in "${MACHINES[@]}"; do MACHINE_LOAD["$m"]=0; done

ACTIVE_MACHINE=(); ACTIVE_PID=(); ACTIVE_REP=(); ACTIVE_SUBFOLDER=(); ACTIVE_LOGFILE=()

# Dispatches the next queued job (JOB_REPLICAS[$QUEUE_IDX] / JOB_ANALYSES[$QUEUE_IDX] /
# JOB_TYPE[$QUEUE_IDX]) onto $1. Always advances QUEUE_IDX exactly once, regardless
# of outcome, so the queue can never stall. Only appends to the ACTIVE_* arrays
# (i.e. only occupies a machine slot) on a confirmed successful dispatch — dry
# runs, failed indexing, and unconfirmed PIDs are all fire-and-forget and free
# their slot immediately.
dispatch_job() {
    local machine="$1"
    local rep="${JOB_REPLICAS[$QUEUE_IDX]}"
    local entry="${JOB_ANALYSES[$QUEUE_IDX]}"
    local type="${JOB_TYPE[$QUEUE_IDX]}"
    local job_num=$(( QUEUE_IDX + 1 ))
    QUEUE_IDX=$(( QUEUE_IDX + 1 ))

    local subfolder="${entry%%:*}"
    local rest="${entry#*:}"
    local script="${rest%%:*}"
    local extras="${rest#*:}"
    [ "$extras" == "$rest" ] && extras=""   # no extras given

    local host; host=$(ssh_host "$machine")
    local ana_dir idx_dir
    if [ "$type" == "pooled" ]; then
        ana_dir="$POOLED_DIR/$subfolder"
    else
        local rep_dir="$PRD_DIR/$rep"
        idx_dir="$rep_dir/analysis/00_index"
        ana_dir="$rep_dir/analysis/$subfolder"
    fi

    log "=== Job $job_num/$NUM_JOBS — ${type} ${rep}, analysis $subfolder → $machine ==="

    mkdir -p "$ana_dir"
    [ "$type" == "pooled" ] || { mkdir -p "$idx_dir"; cp -f "$SCRIPTS_DIR/analysis/00_index/index_maker.sh" "$idx_dir/"; }
    cp -f "$SCRIPTS_DIR/analysis/$subfolder/$script" "$ana_dir/"
    if [ -n "$extras" ]; then
        IFS=',' read -ra EXTRA_FILES <<< "$extras"
        for ef in "${EXTRA_FILES[@]}"; do
            cp -f "$SCRIPTS_DIR/analysis/$subfolder/$ef" "$ana_dir/"
        done
    fi

    if [ "$DRY_RUN" -eq 1 ]; then
        if [ "$type" == "pooled" ]; then
            log "  [DRY RUN] would ssh $host, then nohup $script (pooled, no per-job indexing)"
        else
            log "  [DRY RUN] would ssh $host, build index if needed, then nohup $script"
        fi
        echo "$rep|$subfolder|$machine|$ana_dir|DRYRUN|-" >> "$MANIFEST"
        return
    fi

    # ── INDEX (foreground, blocking, skipped if already built) ──────────────
    # Per-replica jobs only — pooled jobs reuse rep1's index (built the first
    # time any per-replica analysis ran) and don't index themselves.
    # Always uses the default binary — index building (make_ndx) isn't
    # version-pinned the way do_dssp is.
    if [ "$type" != "pooled" ]; then
        if [ "$FORCE_REINDEX" -eq 1 ] || [ ! -f "$idx_dir/index_analysis.ndx" ]; then
            log "  Building index on $machine (foreground)..."
            if ! remote_exec "$host" "cd '$idx_dir' && GROM='$GROM_BIN_DEFAULT' bash index_maker.sh"; then
                log "  ERROR: index_maker.sh failed on $machine for $rep — skipping this job."
                echo "$rep|$subfolder|$machine|$ana_dir|INDEX_FAILED|-" >> "$MANIFEST"
                return
            fi
        else
            log "  Index already present for $rep — skipping (FORCE_REINDEX=0)."
        fi
    fi

    # Resolve which gmx binary this specific analysis needs — falls back to
    # the pool default unless GROM_BIN_OVERRIDE has an entry for $subfolder.
    # Harmless for pooled Python analyses that don't call gmx directly.
    local grom_bin="${GROM_BIN_OVERRIDE[$subfolder]:-$GROM_BIN_DEFAULT}"
    [ -n "${GROM_BIN_OVERRIDE[$subfolder]:-}" ] && log "  Using pinned gmx binary for $subfolder: $grom_bin"

    # ── ANALYSIS (background, nohup, detached) ───────────────────────────────
    local logfile="$ana_dir/${subfolder}_$(date +%Y%m%d_%H%M%S).log"
    log "  Dispatching $script on $machine (nohup, background)..."
    local pid
    local job_cmd="cd '$ana_dir' && GROM='$grom_bin' nohup bash '$script' > '$logfile' 2>&1 < /dev/null"
    pid=$(remote_bg "$host" "$job_cmd")

    if [ -z "$pid" ]; then
        log "  WARNING: could not confirm PID for $rep/$subfolder on $machine — treating as fire-and-forget (not tracked for completion; check manually)."
        echo "$rep|$subfolder|$machine|$ana_dir|UNKNOWN|$logfile" >> "$MANIFEST"
        return
    fi

    log "  → PID $pid on $machine. Log: $logfile"
    echo "$rep|$subfolder|$machine|$ana_dir|$pid|$logfile" >> "$MANIFEST"

    ACTIVE_MACHINE+=("$machine")
    ACTIVE_PID+=("$pid")
    ACTIVE_REP+=("$rep")
    ACTIVE_SUBFOLDER+=("$subfolder")
    ACTIVE_LOGFILE+=("$logfile")
    MACHINE_LOAD["$machine"]=$(( MACHINE_LOAD["$machine"] + 1 ))
}

while [ "$QUEUE_IDX" -lt "$NUM_JOBS" ] || [ "${#ACTIVE_PID[@]}" -gt 0 ]; do
    # Fill every open slot before waiting on anything
    for machine in "${MACHINES[@]}"; do
        while [ "${MACHINE_LOAD[$machine]}" -lt "$MAX_JOBS_PER_MACHINE" ] && [ "$QUEUE_IDX" -lt "$NUM_JOBS" ]; do
            dispatch_job "$machine"
        done
    done

    if [ "$QUEUE_IDX" -ge "$NUM_JOBS" ] && [ "${#ACTIVE_PID[@]}" -eq 0 ]; then
        break   # everything dispatched and confirmed finished/fire-and-forget
    fi

    log "  ${#ACTIVE_PID[@]} job(s) running, $(( NUM_JOBS - QUEUE_IDX )) still queued — checking again in ${POLL_INTERVAL}s..."
    sleep "$POLL_INTERVAL"

    # Poll active jobs; iterate backwards so unset-by-index is safe
    for ((idx=${#ACTIVE_PID[@]}-1; idx>=0; idx--)); do
        machine="${ACTIVE_MACHINE[$idx]}"
        pid="${ACTIVE_PID[$idx]}"
        host=$(ssh_host "$machine")
        if ! remote_alive "$host" "$pid"; then
            log "  ✓ ${ACTIVE_REP[$idx]}/${ACTIVE_SUBFOLDER[$idx]} finished on $machine (PID $pid). Log: ${ACTIVE_LOGFILE[$idx]}"
            MACHINE_LOAD["$machine"]=$(( MACHINE_LOAD["$machine"] - 1 ))
            unset 'ACTIVE_MACHINE[idx]' 'ACTIVE_PID[idx]' 'ACTIVE_REP[idx]' 'ACTIVE_SUBFOLDER[idx]' 'ACTIVE_LOGFILE[idx]'
        fi
    done
    # re-pack arrays to close gaps left by unset
    ACTIVE_MACHINE=("${ACTIVE_MACHINE[@]}")
    ACTIVE_PID=("${ACTIVE_PID[@]}")
    ACTIVE_REP=("${ACTIVE_REP[@]}")
    ACTIVE_SUBFOLDER=("${ACTIVE_SUBFOLDER[@]}")
    ACTIVE_LOGFILE=("${ACTIVE_LOGFILE[@]}")
done

echo
log "✅ All $NUM_JOBS job(s) dispatched and confirmed finished (or fire-and-forget). Manifest: $MANIFEST"
if [ "$EXECUTION_MODE" == "local" ]; then
    log "   Check progress : tail -f <logfile>"
    log "   Check alive    : ps -p <pid>"
else
    log "   Check progress : ssh <machine> tail -f <logfile>"
    log "   Check alive    : ssh <machine> ps -p <pid>"
fi