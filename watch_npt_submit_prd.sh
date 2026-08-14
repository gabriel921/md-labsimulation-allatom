#!/bin/bash
# =============================================================================
# watch_npt_submit_prd.sh  (multi-replicate version)
#
# Monitors the (single, shared) NPT equilibration folder for completion of
# its last segment, then submits production MD for EVERY replicate PRD
# folder given on the command line. One watcher, one NPT wait, N fan-out
# submissions. Each replicate is checked off as soon as it submits
# successfully; if a submission fails, only that replicate is retried on
# the next poll cycle — the others are not touched again.
#
# Usage (standalone):
#   bash watch_npt_submit_prd.sh  NPT_DIR  PRD_DIR_1  [PRD_DIR_2 ... PRD_DIR_N]
#
# Backward compatible: called with exactly one PRD dir, it behaves exactly
# like the single-replicate version.
#
# Launched from the master setup script as:
#   nohup bash watch_npt_submit_prd.sh "$NPT" "${PRD_DIRS[@]}" &
#   WATCHER_PID=$!
#
# Resumability: each replicate dir gets a ".prd_submitted" marker once its
# job is submitted, so re-running the watcher (e.g. after a crash) will not
# resubmit replicates that already went out.
# =============================================================================

# ── ARGUMENTS ─────────────────────────────────────────────────────────────
NPT_DIR="$1"
shift || true
PRD_DIRS=("$@")

POLL_INTERVAL=20   # seconds between checks (NPT wait, and retry cycles)
SETTLE_WAIT=5       # seconds to wait after the file appears before submitting
MAX_RETRIES=10      # per-replicate retry attempts if submission fails

# Filename used to stage the NPT endpoint structure inside each replica dir.
# MUST NOT look like a production segment name (001.gro, 002.gro, ...) —
# slurm_MD_prd_allatom.sh detects "next segment to run" by checking for the
# existence of 001.gro, 002.gro, etc. in the replica directory. If the NPT
# endpoint happened to be copied in under one of those exact names (this
# happens whenever NPT's own last segment is 001, i.e. NPT_SEGMENTS=1), the
# production script mistakes the copied starting structure for an already
# -finished segment and silently starts numbering its real output one
# segment later than intended (e.g. first real output is 002.xtc, not
# 001.xtc) — which also throws off the Segments= self-requeue count by one.
NPT_ENDPOINT_NAME="npt_final.gro"

# ── HELPERS ──────────────────────────────────────────────────────────────────
log() { echo -e "[$(date '+%H:%M:%S')] [WATCHER] $*"; }
die() { echo -e "[$(date '+%H:%M:%S')] [WATCHER] ERROR: $*" >&2; exit 1; }

# ── VALIDATE ARGUMENTS ────────────────────────────────────────────────────────
[ -n "$NPT_DIR" ] || die "No NPT directory supplied. Usage: $0 NPT_DIR PRD_DIR_1 [PRD_DIR_2 ...]"
[ "${#PRD_DIRS[@]}" -ge 1 ] || die "No PRD directories supplied. Usage: $0 NPT_DIR PRD_DIR_1 [PRD_DIR_2 ...]"
[ -d "$NPT_DIR" ] || die "NPT directory not found: $NPT_DIR"

for d in "${PRD_DIRS[@]}"; do
    [ -d "$d" ]                                  || die "PRD directory not found: $d"
    [ -f "$d/slurm_MD_prd_allatom.dat" ]         || die "PRD dat file not found: $d/slurm_MD_prd_allatom.dat"
    [ -f "$d/slurm_MD_prd_allatom.sh"  ]         || die "PRD submit script not found: $d/slurm_MD_prd_allatom.sh"
done

NPT_DAT="$NPT_DIR/slurm_MD_eq_NPT_allatom.dat"
[ -f "$NPT_DAT" ] || die "NPT dat file not found: $NPT_DAT"

# ── READ SEGMENT COUNT FROM NPT DAT ──────────────────────────────────────────
NPT_SEGMENTS=$(grep -E '^[[:space:]]*Segments=' "$NPT_DAT" | tail -1 | cut -d= -f2 | tr -d '[:space:]')
[ -z "$NPT_SEGMENTS" ] && die "Could not parse Segments= from $NPT_DAT"

LAST_SEG=$(printf "%03d" "$NPT_SEGMENTS")
TARGET_XTC="$NPT_DIR/${LAST_SEG}.xtc"
TARGET_GRO="$NPT_DIR/${LAST_SEG}.gro"

log "NPT directory     : $NPT_DIR"
log "Replicate PRD dirs: ${PRD_DIRS[*]}"
log "Waiting for       : $TARGET_XTC  (segment $LAST_SEG of $NPT_SEGMENTS)"
log "Poll interval     : ${POLL_INTERVAL}s   Settle wait: ${SETTLE_WAIT}s"
log "----------------------------------------------------------------------"

# ── PHASE 1 — WAIT FOR NPT TO FINISH (single shared wait) ────────────────────
while [ ! -f "$TARGET_XTC" ]; do
    log "NPT not finished yet — next check in ${POLL_INTERVAL}s..."
    sleep "$POLL_INTERVAL"
done

log "Found $TARGET_XTC — waiting ${SETTLE_WAIT}s for I/O to settle..."
sleep "$SETTLE_WAIT"

if [ ! -f "$TARGET_GRO" ]; then
    log "WARNING: $TARGET_GRO not yet present — waiting one more poll cycle..."
    sleep "$POLL_INTERVAL"
    [ -f "$TARGET_GRO" ] || die "$TARGET_GRO still missing after extra wait. Aborting."
fi
log "$TARGET_GRO confirmed present. NPT is done — fanning out to ${#PRD_DIRS[@]} replicate(s)."
log "----------------------------------------------------------------------"

# ── PHASE 2 — SUBMIT PRODUCTION FOR EACH REPLICATE ───────────────────────────
# Track retry counts per replicate dir.
declare -A RETRIES
for d in "${PRD_DIRS[@]}"; do RETRIES["$d"]=0; done

# Working list of replicates still pending. Anything already marked
# .prd_submitted (e.g. from a previous watcher run) is skipped up front.
PENDING=()
for d in "${PRD_DIRS[@]}"; do
    if [ -f "$d/.prd_submitted" ]; then
        log "Replica already submitted previously, skipping: $d"
    else
        PENDING+=("$d")
    fi
done

submit_one() {
    local d="$1"

    # Stage the NPT endpoint under a name that can NEVER collide with
    # production's own segment numbering (001, 002, ...) — see
    # NPT_ENDPOINT_NAME comment above.
    cp -f "$TARGET_GRO" "$d/${NPT_ENDPOINT_NAME}" || { log "  copy of NPT endpoint failed for $d"; return 1; }

    sed -i "s|^[[:space:]]*gro=.*|gro=${NPT_ENDPOINT_NAME}|" "$d/slurm_MD_prd_allatom.dat" \
        || { log "  failed to patch gro= in $d/slurm_MD_prd_allatom.dat"; return 1; }

    log "Submitting production for replica: $d"
    ( cd "$d" && bash slurm_MD_prd_allatom.sh ) || { log "  sbatch submission failed for $d"; return 1; }

    touch "$d/.prd_submitted"
    log "  Replica objective reached (job submitted): $d — removing from watch list."
    return 0
}

while [ "${#PENDING[@]}" -gt 0 ]; do
    STILL_PENDING=()
    for d in "${PENDING[@]}"; do
        if submit_one "$d"; then
            continue
        fi

        RETRIES["$d"]=$(( RETRIES["$d"] + 1 ))
        if [ "${RETRIES[$d]}" -ge "$MAX_RETRIES" ]; then
            log "  Giving up on $d after ${MAX_RETRIES} failed attempts — check it manually."
        else
            log "  Will retry $d (attempt ${RETRIES[$d]}/${MAX_RETRIES}) next cycle."
            STILL_PENDING+=("$d")
        fi
    done

    PENDING=("${STILL_PENDING[@]}")
    if [ "${#PENDING[@]}" -gt 0 ]; then
        log "Remaining replicates to submit: ${PENDING[*]}"
        sleep "$POLL_INTERVAL"
    fi
done

log "----------------------------------------------------------------------"
log "Phase 2 complete — all replicates submitted (or exhausted retries)."
log "----------------------------------------------------------------------"

# ── PHASE 3 — WAIT FOR EACH REPLICA'S PRODUCTION TO FINISH, THEN CONCATENATE ──
# Same fan-out/check-off pattern as Phase 2: one shared watcher, each
# replica dropped from the list once its trajectory is concatenated.
# "Finished" is detected the same way as NPT was — the last segment's .gro
# appearing in that replica's own folder.

CONCAT_SCRIPT_NAME="concat_trj_replica.sh"

declare -A CONCAT_RETRIES
for d in "${PRD_DIRS[@]}"; do CONCAT_RETRIES["$d"]=0; done

CONCAT_PENDING=()
for d in "${PRD_DIRS[@]}"; do
    if [ -f "$d/.concat_done" ]; then
        log "Replica already concatenated previously, skipping: $d"
    elif [ ! -f "$d/.prd_submitted" ]; then
        log "WARNING: production was never submitted for $d — skipping concat wait for it."
    else
        CONCAT_PENDING+=("$d")
    fi
done

last_segment_for() {
    local d="$1"
    local dat="$d/slurm_MD_prd_allatom.dat"
    [ -f "$dat" ] || { echo ""; return; }
    local segs
    segs=$(grep -E '^[[:space:]]*Segments=' "$dat" | tail -1 | cut -d= -f2 | tr -d '[:space:]')
    [ -z "$segs" ] && { echo ""; return; }
    printf "%03d" "$segs"
}

concat_one() {
    local d="$1"
    local last_seg
    last_seg=$(last_segment_for "$d")
    if [ -z "$last_seg" ]; then
        log "  Could not determine segment count for $d (missing/unparsable .dat) — will retry."
        return 1
    fi

    if [ ! -f "$d/${last_seg}.gro" ]; then
        return 2   # production not finished yet — not a failure, just not ready
    fi

    log "Production finished for $d (segment $last_seg) — concatenating..."
    sleep "$SETTLE_WAIT"

    if [ ! -f "$d/$CONCAT_SCRIPT_NAME" ]; then
        log "  ERROR: $CONCAT_SCRIPT_NAME not found in $d — copy it there first."
        return 1
    fi

    ( cd "$d" && bash "$CONCAT_SCRIPT_NAME" ) || { log "  concat script failed for $d"; return 1; }

    touch "$d/.concat_done"
    log "  Replica concatenation objective reached: $d — removing from watch list."
    return 0
}

while [ "${#CONCAT_PENDING[@]}" -gt 0 ]; do
    STILL_PENDING=()
    for d in "${CONCAT_PENDING[@]}"; do
        concat_one "$d"
        status=$?

        if [ "$status" -eq 0 ]; then
            continue
        elif [ "$status" -eq 2 ]; then
            STILL_PENDING+=("$d")   # not finished yet — no retry penalty
            continue
        fi

        CONCAT_RETRIES["$d"]=$(( CONCAT_RETRIES["$d"] + 1 ))
        if [ "${CONCAT_RETRIES[$d]}" -ge "$MAX_RETRIES" ]; then
            log "  Giving up concatenation for $d after ${MAX_RETRIES} failed attempts — check it manually."
        else
            log "  Will retry concatenation for $d (attempt ${CONCAT_RETRIES[$d]}/${MAX_RETRIES}) next cycle."
            STILL_PENDING+=("$d")
        fi
    done

    CONCAT_PENDING=("${STILL_PENDING[@]}")
    if [ "${#CONCAT_PENDING[@]}" -gt 0 ]; then
        sleep "$POLL_INTERVAL"
    fi
done

log "----------------------------------------------------------------------"
log "All replicates concatenated (or exhausted retries). Watcher exiting."
exit 0