#!/bin/bash -e
# =============================================================================
# run_setup.sh
# Master setup script for Aβ42 all-atom MD simulation (CHARMM36m)
#
# Runs sequentially:
#   0) Topology generation   (pdb2gmx)
#   1) Box building          (editconf)
#   2) Solvation             (solvate)
#   3) Ionisation            (genion)
#   4) Energy minimisation   (min_allatom.sh)
#   5) NVT equilibration     (eq_NVT_allatom.sh)
#   6) NPT equilibration     (SLURM by default, or local sequential — see
#   7) Production MD          EXECUTION_MODE below)
#
# Set EXECUTION_MODE=slurm (default) or EXECUTION_MODE=local in the USER
# CONFIGURATION section to choose how steps 6-7 run. Steps 0-5 are identical
# either way.
# =============================================================================
# ── LOG EVERYTHING TO FILE ───────────────────────────────────────────────────
#LOGFILE="run_setup_$(date +%Y%m%d_%H%M%S).log"
#exec >> "$LOGFILE" 2>&1
#echo "============================================================"
#echo " Logging started: $(date)"
#echo " Log file: $LOGFILE"
#echo "============================================================"

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
PDB=ab42_1z0q.pdb                        # Input PDB file (must be in the same directory)
SYS=ab42                                 # System name prefix for output files
GMX=/gromacs/gromacs-2021.2/bin/gmx
export GMXLIB=$HOME/Desktop/GMXLIB       # GROMACS force field library path

BOX_PADDING=1.2                          # Minimum distance from solute to box edge (nm)
BOX_TYPE=cubic                           # Box geometry
ION_CONC=0.1                             # NaCl concentration (M)
NTHREADS=8                               # CPU threads for local steps (min, NVT)
MAXWARN=1                                # Allowed grompp warnings

MDP_DIR=$(pwd)/mdp                       # Directory containing all .mdp files
SCRIPTS_DIR=$(pwd)                       # System root — watcher script lives here
INP_DIR=$(pwd)/input/simulation          # Directory containing .dat files + slurm_*.sh + concat script

FF=13                                    # 13 = CHARMM36m July 2021. Change if needed — check your GMXLIB
WM=1

NUM_CHAINS=1                             #Number of chains in the PDB
TERMINI_N=0                              #Charged = 0, Neutral = 1 (repeats for every chain)
TERMINI_C=0                              #Pay attention to special cases. Check the README file

POSRES_NPT=yes                           # POSRES_NPT: position restraints during NPT equilibration. 
                                         # yes = keep "define = -DPOSRES" in the NPT .mdp     no = comment it out
REPLICATES=3                             # Number of independent production replicates
                                         # Each gets its own 4_prd/repN subfolder,
                                         # all starting from the same NPT endpoint.

# ── EXECUTION MODE ────────────────────────────────────────────────────────────
# "slurm" (default) — unchanged original behaviour: NPT and production run as
#                      self-requeuing SLURM jobs (steps 6-7 below), watched in
#                      the background by watch_npt_submit_prd.sh.
# "local"            — no SLURM involved. NPT then each production replicate
#                       run sequentially, in the foreground, in this shell,
#                       reusing the same grompp/mdrun commands and the same
#                       segment/continuation logic (see run_segments_local.sh).
EXECUTION_MODE=slurm

# ── NPT / PRODUCTION STEP-COUNT AND BLOCK-COUNT CONFIGURATION ────────────────
# These control the `nsteps`/`dt` written into the copies of eq_npt_allatom.mdp
# and prd_allatom.mdp used for the run, and the number of NPT/production
# blocks (SLURM "Segments", or local sequential iterations). Defaults below
# match the values shipped in the original .mdp/.dat files, so leaving these
# untouched reproduces the original behaviour exactly, in either execution mode.
NPT_NSTEPS=50                            # nsteps written into the NPT .mdp for this run
NPT_DT=0.002                             # dt (ps) written into the NPT .mdp for this run
NPT_SEGMENTS=1                           # number of NPT blocks (SLURM Segments= / local loop count)

PROD_NSTEPS=500                          # nsteps written into the production .mdp for this run
PROD_DT=0.002                            # dt (ps) written into the production .mdp for this run
PROD_SEGMENTS=3                          # number of production blocks per replicate

# ── ANALYSIS INDEX GROUPS (system_groups.cfg) ────────────────────────────────
# Written once to $PREP/system_groups.cfg and sourced by every analysis
# script (index_maker.sh, Rg, Sasa, etc.) — this is the ONLY place chain and
# peptide atom ranges should need editing when you move to a different
# protein/system. See $PREP/system_groups.cfg after setup for the full
# format and comments.
ANALYSIS_BACKBONE_GROUP_NAME="Backbone"
ANALYSIS_CHAIN_DEFS=("ChainA:1-627")      # 1-627 is the number of atoms of Aβ42, change accordingly to other systems
ANALYSIS_PEPTIDE_DEFS=()
ANALYSIS_BACKBONE_FOR=(ChainA)
ANALYSIS_ALL_CHAINS_NAME="All_Chains"

ANALYSIS_CONCAT_HAS_SOLUTE=1              # 1 = center on a solute; 0 = PBC only
ANALYSIS_CONCAT_SOLUTE_GROUP="Protein"    # group to center on (if above = 1)
ANALYSIS_CONCAT_CENTER_OUTPUT_GROUP="Protein"  # group written to the final trajectory
# ── END USER CONFIGURATION ───────────────────────────────────────────────────

# ── FOLDER STRUCTURE ─────────────────────────────────────────────────────────
PREP=0_prep
MIN=1_min
NVT=2_nvt
NPT=3_npt
PRD=4_prd

TOP=topol.top
POS=posre.itp
NDX=index.ndx

# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────
log()        { echo -e "\n[$(date '+%H:%M:%S')] >>> $*"; }
die()        { echo -e "\n[ERROR] $*" >&2; exit 1; }
check_file() { [ -f "$1" ] || die "Required file not found: $1"; }
check_cmd()  { command -v "$1" &>/dev/null || die "Command not found: $1 — check your GROMACS installation."; }

# Patches nsteps/dt into a copied .mdp file (never touches the template in
# $MDP_DIR — only the per-run copy passed as $1).
patch_mdp_steps() {
    local mdp_file="$1" nsteps="$2" dt="$3"
    sed -i "s/^nsteps[[:space:]]*=.*/nsteps      = ${nsteps}          ; set by run_setup.sh/" "$mdp_file"
    sed -i "s/^dt[[:space:]]*=.*/dt          = ${dt}             ; set by run_setup.sh/"     "$mdp_file"
}

# <<< CHANGED (new function) — comments out the define line in the per-run copy only
patch_mdp_posres() {
    local mdp_file="$1" use_posres="$2"
    if [ "$use_posres" = "no" ]; then
        sed -i 's/^define[[:space:]]*=.*/; define disabled by run_setup.sh (POSRES_NPT=no)/' "$mdp_file"
    fi
}

# Patches Segments= into a copied SLURM .dat file (never touches the
# template in $INP_DIR — only the per-run copy passed as $1).
patch_dat_segments() {
    local dat_file="$1" segments="$2"
    sed -i "s/^Segments=.*/Segments=${segments}/" "$dat_file"
}

# Checks that a requested GMXLIB directory actually exists (a clear, common
# misconfiguration), and — best-effort only — lists the .ff directories it
# and the GROMACS installation can see, as a sanity aid. Only the missing-
# GMXLIB-directory case is fatal; pdb2gmx itself remains the final authority
# on whether FF=<n> is actually valid, so this check is not allowed to block
# a run that pdb2gmx would otherwise have handled fine (e.g. non-standard
# GROMACS install layouts where the share/top path guess below doesn't hold).
check_forcefield() {
    local ff_index="$1"
    local gmx_share
    gmx_share="$(dirname "$(dirname "$GMX")")/share/gromacs/top"

    if [ -n "$GMXLIB" ] && [ ! -d "$GMXLIB" ]; then
        die "GMXLIB is set to '$GMXLIB' but that directory does not exist. Fix GMXLIB in run_setup.sh, or unset it to use only the force fields bundled with GROMACS."
    fi

    local -a ff_list=()
    if [ -n "$GMXLIB" ] && [ -d "$GMXLIB" ]; then
        while IFS= read -r -d '' d; do ff_list+=("$(basename "$d")"); done \
            < <(find "$GMXLIB" -maxdepth 1 -type d -name '*.ff' -print0 2>/dev/null | sort -z)
    fi
    if [ -d "$gmx_share" ]; then
        while IFS= read -r -d '' d; do ff_list+=("$(basename "$d")"); done \
            < <(find "$gmx_share" -maxdepth 1 -type d -name '*.ff' -print0 2>/dev/null | sort -z)
    fi

    if [ "${#ff_list[@]}" -eq 0 ]; then
        log "WARNING: could not list any force fields under GMXLIB ('$GMXLIB') or the guessed GROMACS path ('$gmx_share'). Continuing anyway — pdb2gmx will report a clear error at STEP 0 if FF=$ff_index is genuinely invalid."
        return 0
    fi

    if ! [[ "$ff_index" =~ ^[0-9]+$ ]] || [ "$ff_index" -lt 1 ] || [ "$ff_index" -gt "${#ff_list[@]}" ]; then
        {
            echo "WARNING: FF=$ff_index does not match any of the force fields this script could find."
            echo "Available force fields (GMXLIB first, then the GROMACS installation) — for reference only, numbering here may not exactly match pdb2gmx's own menu:"
            local n=1
            for f in "${ff_list[@]}"; do echo "    $n) $f"; n=$((n + 1)); done
            echo "If STEP 0 (pdb2gmx) fails, set FF to a valid menu index and re-run."
        } >&2
        return 0
    fi
    log "Force field check OK — FF=$ff_index -> ${ff_list[$((ff_index - 1))]}"
}

# ── PREFLIGHT CHECKS ─────────────────────────────────────────────────────────
log "Running preflight checks..."
check_cmd  "$GMX"
check_file "$PDB"
check_file "$MDP_DIR/ions_allatom.mdp"
check_file "$MDP_DIR/min_allatom.mdp"
check_file "$MDP_DIR/eq_nvt_allatom.mdp"
check_file "$MDP_DIR/eq_npt_allatom.mdp"
check_file "$MDP_DIR/prd_allatom.mdp"
check_file "$INP_DIR/min_allatom.sh"
check_file "$INP_DIR/eq_NVT_allatom.sh"
check_file "$INP_DIR/concat_trj_replica.sh"

case "$EXECUTION_MODE" in
    slurm)
        check_file "$INP_DIR/slurm_MD_eq_NPT_allatom.sh"
        check_file "$INP_DIR/slurm_MD_eq_NPT_allatom.dat"
        check_file "$INP_DIR/slurm_MD_prd_allatom.sh"
        check_file "$INP_DIR/slurm_MD_prd_allatom.dat"
        check_file "$SCRIPTS_DIR/watch_npt_submit_prd.sh"
        ;;
    local)
        check_file "$INP_DIR/run_segments_local.sh"
        ;;
    *)
        die "Unknown EXECUTION_MODE='$EXECUTION_MODE' — must be 'slurm' or 'local'."
        ;;
esac

# <<< CHANGED (new block) — validate the two new options
[[ "$NUM_CHAINS" =~ ^[1-9][0-9]*$ ]] || die "NUM_CHAINS='$NUM_CHAINS' — must be a positive integer."
for t in "$TERMINI_N" "$TERMINI_C"; do
    [[ "$t" =~ ^[0-9]+$ ]] || die "Invalid TERMINI value '$t' — must be a pdb2gmx menu number (0, 1, 2, ...)."
done
case "$POSRES_NPT" in
    yes|no) ;;
    *) die "Unknown POSRES_NPT='$POSRES_NPT' — must be 'yes' or 'no'." ;;
esac

check_forcefield "$FF"
log "All checks passed."

# ── CREATE DIRECTORY STRUCTURE ───────────────────────────────────────────────
log "Creating directory structure..."
mkdir -p $PREP $MIN $NVT $NPT $PRD
echo "  $PREP/  — system preparation (topology, box, solvation, ions)"
echo "  $MIN/   — energy minimisation"
echo "  $NVT/   — NVT equilibration"
echo "  $NPT/   — NPT equilibration ($EXECUTION_MODE mode)"
echo "  $PRD/   — production MD ($EXECUTION_MODE mode)"

# ── WRITE system_groups.cfg ───────────────────────────────────────────────────
# Single source of truth for chain/peptide atom ranges, consumed by
# index_maker.sh and any later analysis script. Regenerated every run from
# the ANALYSIS_* arrays above, so editing those and re-running setup is all
# that's needed to adapt to a different protein.
log "Writing $PREP/system_groups.cfg..."
{
    echo "# Auto-generated by run_setup.sh on $(date) — edit the ANALYSIS_*"
    echo "# arrays at the top of that script and re-run setup to regenerate."
    echo
    echo "BACKBONE_GROUP_NAME=\"${ANALYSIS_BACKBONE_GROUP_NAME}\""
    echo
    echo "CHAIN_DEFS=("
    for def in "${ANALYSIS_CHAIN_DEFS[@]}"; do echo "    \"${def}\""; done
    echo ")"
    echo
    echo "PEPTIDE_DEFS=("
    for def in "${ANALYSIS_PEPTIDE_DEFS[@]}"; do echo "    \"${def}\""; done
    echo ")"
    echo
    echo "BACKBONE_FOR=(${ANALYSIS_BACKBONE_FOR[*]})"
    echo
    echo "ALL_CHAINS_NAME=\"${ANALYSIS_ALL_CHAINS_NAME}\""
    echo
    echo "CONCAT_HAS_SOLUTE=${ANALYSIS_CONCAT_HAS_SOLUTE}"
    echo "CONCAT_SOLUTE_GROUP=\"${ANALYSIS_CONCAT_SOLUTE_GROUP}\""
    echo "CONCAT_CENTER_OUTPUT_GROUP=\"${ANALYSIS_CONCAT_CENTER_OUTPUT_GROUP}\""
} > "$PREP/system_groups.cfg"
log "$PREP/system_groups.cfg written."

# ── STEP 0 — Topology generation ─────────────────────────────────────────────
log "STEP 0 — pdb2gmx: generating topology and .gro from $PDB..."
cd $PREP
{ printf "%s\n%s\n" "$FF" "$WM"
  for i in $(seq 1 "$NUM_CHAINS"); do printf "%s\n%s\n" "$TERMINI_N" "$TERMINI_C"; done; } |   # <<< CHANGED (N/C answers repeated per chain)
    $GMX pdb2gmx -f ../$PDB -ignh -ter -o ${SYS}.gro                      # <<< CHANGED (-ter added)
check_file "${SYS}.gro"
check_file "$TOP"

# ── STEP 1 — Build simulation box ────────────────────────────────────────────
log "STEP 1 — editconf: building $BOX_TYPE box (padding = ${BOX_PADDING} nm)..."
$GMX editconf \
    -f  ${SYS}.gro \
    -o  ${SYS}_box.gro \
    -c \
    -d  $BOX_PADDING \
    -bt $BOX_TYPE
check_file "${SYS}_box.gro"

# ── STEP 2 — Solvation ───────────────────────────────────────────────────────
log "STEP 2 — solvate: filling box with water (SPC/E)..."
$GMX solvate \
    -cp ${SYS}_box.gro \
    -p  $TOP \
    -cs spc216.gro \
    -o  ${SYS}_solv.gro
check_file "${SYS}_solv.gro"

# ── STEP 3 — Add ions ────────────────────────────────────────────────────────
log "STEP 3 — genion: neutralising system and adding ${ION_CONC} M NaCl..."
$GMX grompp \
    -f $MDP_DIR/ions_allatom.mdp \
    -c ${SYS}_solv.gro \
    -p $TOP \
    -o ions.tpr \
    -maxwarn $MAXWARN

echo "SOL" | $GMX genion \
    -s ions.tpr \
    -o ${SYS}_solv_ions.gro \
    -p $TOP \
    -pname SOD \
    -nname CLA \
    -neutral \
    -conc $ION_CONC
check_file "${SYS}_solv_ions.gro"
cd ..

# ── STEP 4 — Energy minimisation ─────────────────────────────────────────────
log "STEP 4 — Energy minimisation..."
cd $MIN
GMX=$GMX \
NTHREADS=$NTHREADS \
MDP=$MDP_DIR/min_allatom.mdp \
GRO=../$PREP/${SYS}_solv_ions.gro \
TOP=../$PREP/$TOP \
    bash $INP_DIR/min_allatom.sh
check_file "em.gro"
cd ..

# ── STEP 5 — NVT equilibration ───────────────────────────────────────────────
log "STEP 5 — NVT equilibration..."
cd $NVT
GMX=$GMX \
NTHREADS=$NTHREADS \
MDP=$MDP_DIR/eq_nvt_allatom.mdp \
GRO=../$MIN/em.gro \
TOP=../$PREP/$TOP \
NDX=$NDX \
    bash $INP_DIR/eq_NVT_allatom.sh
check_file "nvt.gro"
cd ..

if [ "$EXECUTION_MODE" = "slurm" ]; then
    # ── STEP 6 — NPT equilibration (SLURM) ───────────────────────────────────
    log "STEP 6 — Setting up NPT folder and submitting SLURM job..."
    cd $NPT
    cp ../$PREP/$TOP       .
    cp ../$PREP/$POS       .
    cp ../$NVT/$NDX        .
    cp ../$NVT/nvt.gro     .
    cp ../$NVT/nvt.cpt     .
    cp $MDP_DIR/eq_npt_allatom.mdp .
    patch_mdp_steps eq_npt_allatom.mdp "$NPT_NSTEPS" "$NPT_DT"
    patch_mdp_posres eq_npt_allatom.mdp "$POSRES_NPT"                     # <<< CHANGED (new line)
    cp $INP_DIR/slurm_MD_eq_NPT_allatom.dat .
    patch_dat_segments slurm_MD_eq_NPT_allatom.dat "$NPT_SEGMENTS"
    cp $INP_DIR/slurm_MD_eq_NPT_allatom.sh  .
    chmod +x slurm_MD_eq_NPT_allatom.sh
    ./slurm_MD_eq_NPT_allatom.sh
    log "NPT SLURM job submitted from $(pwd)."
    cd ..

    # ── STEP 7 — Production MD folders (one per replicate) + background watcher ──
    log "STEP 7 — Preparing $REPLICATES production replicate folder(s)..."
    cd $PRD

    PRD_REP_DIRS=()   # absolute paths, filled in below, passed to the watcher

    for r in $(seq -w 1 "$REPLICATES"); do
        REPDIR="rep${r}"
        mkdir -p "$REPDIR"
        cp ../$PREP/$TOP               "$REPDIR/"
        cp ../$NPT/$NDX                "$REPDIR/"
        cp $MDP_DIR/prd_allatom.mdp           "$REPDIR/"
        patch_mdp_steps "$REPDIR/prd_allatom.mdp" "$PROD_NSTEPS" "$PROD_DT"
        cp $INP_DIR/slurm_MD_prd_allatom.dat  "$REPDIR/"
        patch_dat_segments "$REPDIR/slurm_MD_prd_allatom.dat" "$PROD_SEGMENTS"
        cp $INP_DIR/slurm_MD_prd_allatom.sh   "$REPDIR/"
        cp $INP_DIR/concat_trj_replica.sh     "$REPDIR/"

        # Suffix the job Name so each replicate shows up distinctly in squeue,
        # e.g. Name=ab42_prd  ->  Name=ab42_prd_rep1
        if grep -qE '^[[:space:]]*Name=' "$REPDIR/slurm_MD_prd_allatom.dat"; then
            sed -i -E "s|^([[:space:]]*Name=[^[:space:]#]+)|\1_rep${r}|" "$REPDIR/slurm_MD_prd_allatom.dat"
        else
            log "WARNING: no Name= line found in $REPDIR/slurm_MD_prd_allatom.dat — job names will collide in squeue."
        fi

        chmod +x "$REPDIR/slurm_MD_prd_allatom.sh"
        PRD_REP_DIRS+=("$(pwd)/$REPDIR")
        echo "  $PRD/$REPDIR/  — replicate $r"
    done
    log "$REPLICATES production folder(s) ready."
    cd ..

    # Launch a single watcher in the background, pointed at the absolute NPT path
    # and every replicate PRD path, so it survives any later cd changes in the
    # calling shell. It waits once for NPT to finish, then fans out production
    # submissions across all replicates, checking each off as it submits.
    WATCHER_LOG="watch_npt_submit_prd_$(date +%Y%m%d_%H%M%S).log"
    nohup bash "$SCRIPTS_DIR/watch_npt_submit_prd.sh" \
          "$(pwd)/$NPT" \
          "${PRD_REP_DIRS[@]}" \
          >> "$WATCHER_LOG" 2>&1 &
    WATCHER_PID=$!

    # ── SUMMARY ──────────────────────────────────────────────────────────────────
    echo ""
    echo "============================================================"
    echo "  Setup complete — $(date)"
    echo "============================================================"
    echo ""
    echo "  Directory layout:"
    echo "    $PREP/  — topology, solvated+ionised structure"
    echo "    $MIN/   — energy minimisation output (em.gro)"
    echo "    $NVT/   — NVT equilibration output  (nvt.gro)"
    echo "    $NPT/   — NPT equilibration         (SLURM job submitted, $NPT_SEGMENTS segment(s))"
    echo "    $PRD/   — production MD             ($REPLICATES replicate(s), auto-submitted, $PROD_SEGMENTS segment(s) each)"
    for r in $(seq -w 1 "$REPLICATES"); do
        echo "      $PRD/rep${r}/"
    done
    echo ""
    echo "  ✔  Background watcher started (PID $WATCHER_PID)"
    echo "     It will poll $NPT/ every 20 s for the last NPT segment .xtc,"
    echo "     then submit production for all $REPLICATES replicate(s):"
    for r in $(seq -w 1 "$REPLICATES"); do
        echo "       - $PRD/rep${r}/slurm_MD_prd_allatom.sh"
    done
    echo "     Once each replicate's production finishes, the same watcher"
    echo "     automatically concatenates its trajectory (concat_trj_replica.sh)"
    echo "     into $PRD/repN/traj_all_center_PBC.xtc — the file the analysis"
    echo "     scripts (Rg, etc.) expect to find."
    echo "     Watcher log: $WATCHER_LOG"
    echo ""
    echo "  To monitor the watcher:"
    echo "     tail -f $WATCHER_LOG"
    echo ""
    echo "  To cancel the watcher (if needed):"
    echo "     kill $WATCHER_PID"
    echo ""

else
    # ── LOCAL MODE — STEP 6 — NPT equilibration (sequential, foreground) ────────
    log "STEP 6 — NPT equilibration (LOCAL mode, $NPT_SEGMENTS segment(s))..."
    cd $NPT
    cp ../$PREP/$TOP       .
    cp ../$NVT/$NDX        .
    cp ../$NVT/nvt.gro     .
    cp ../$NVT/nvt.cpt     .
    cp $MDP_DIR/eq_npt_allatom.mdp .
    patch_mdp_steps eq_npt_allatom.mdp "$NPT_NSTEPS" "$NPT_DT"
    patch_mdp_posres eq_npt_allatom.mdp "$POSRES_NPT"                     # <<< CHANGED (new line)

    GMX=$GMX NTHREADS=$NTHREADS MDP=eq_npt_allatom.mdp TOP=$TOP NDX=$NDX \
    GRO=nvt.gro SEGMENTS=$NPT_SEGMENTS MAXWARN=$MAXWARN \
        bash $INP_DIR/run_segments_local.sh

    NPT_LAST_SEG=$(printf "%03d" "$NPT_SEGMENTS")
    check_file "${NPT_LAST_SEG}.gro"
    NPT_FINAL_GRO="$(pwd)/${NPT_LAST_SEG}.gro"
    log "NPT complete (local): $NPT_FINAL_GRO"
    cd ..

    # ── LOCAL MODE — STEP 7 — Production MD, replicates run one after another ──
    log "STEP 7 — Production MD (LOCAL mode, $REPLICATES replicate(s), $PROD_SEGMENTS segment(s) each)..."
    cd $PRD

    for r in $(seq -w 1 "$REPLICATES"); do
        REPDIR="rep${r}"
        mkdir -p "$REPDIR"
        cp ../$PREP/$TOP        "$REPDIR/"
        cp ../$NPT/$NDX         "$REPDIR/"
        cp $MDP_DIR/prd_allatom.mdp "$REPDIR/"
        patch_mdp_steps "$REPDIR/prd_allatom.mdp" "$PROD_NSTEPS" "$PROD_DT"
        cp $INP_DIR/concat_trj_replica.sh "$REPDIR/"
        cp "$NPT_FINAL_GRO" "$REPDIR/npt_final.gro"

        log "  Replicate $r — running $PROD_SEGMENTS segment(s) sequentially..."
        (
            cd "$REPDIR"
            GMX=$GMX NTHREADS=$NTHREADS MDP=prd_allatom.mdp TOP=$TOP NDX=$NDX \
            GRO=npt_final.gro SEGMENTS=$PROD_SEGMENTS MAXWARN=$MAXWARN \
                bash $INP_DIR/run_segments_local.sh
        )

        log "  Replicate $r — concatenating trajectory..."
        ( cd "$REPDIR" && bash concat_trj_replica.sh )

        echo "  $PRD/$REPDIR/  — replicate $r complete (local)"
    done
    cd ..

    # ── SUMMARY ──────────────────────────────────────────────────────────────
    echo ""
    echo "============================================================"
    echo "  Setup complete — $(date)  [LOCAL execution mode]"
    echo "============================================================"
    echo ""
    echo "  Directory layout:"
    echo "    $PREP/  — topology, solvated+ionised structure"
    echo "    $MIN/   — energy minimisation output (em.gro)"
    echo "    $NVT/   — NVT equilibration output  (nvt.gro)"
    echo "    $NPT/   — NPT equilibration         ($NPT_SEGMENTS segment(s), completed)"
    echo "    $PRD/   — production MD             ($REPLICATES replicate(s), $PROD_SEGMENTS segment(s) each, completed)"
    for r in $(seq -w 1 "$REPLICATES"); do
        echo "      $PRD/rep${r}/traj_all_center_PBC.xtc  — ready for run_setup_analysis.sh"
    done
    echo ""
    echo "  All stages ran sequentially in this shell — no SLURM was used."
    echo ""
fi