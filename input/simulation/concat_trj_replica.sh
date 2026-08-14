#!/bin/bash -e
# =============================================================================
# concat_trj_replica.sh
#
# Concatenates all finished production segments for the replica this script
# is run from into one PBC-corrected (optionally solute-centered) trajectory,
# plus a skip-reduced version for VMD/PyMOL visualization. Run this from
# inside a replica's own production folder, e.g. 4_prd/rep1/, once all
# segments have finished (the watcher does this automatically).
#
# Outputs (fixed names — no duration embedded in the filename, since
# segment count/length varies per run and embedding it made downstream
# analysis scripts fragile about guessing the right file):
#   traj_all_center_PBC.xtc      <- the trajectory analysis scripts use
#   traj_all_center_PBC_VMD.xtc  <- skip-reduced, for visual inspection
# traj_all_no_PBC.xtc (intermediate) is deleted automatically.
#
# Env vars (all optional, fallback defaults match the standard layout):
#   GROM    - gmx binary                    (default: /gromacs/gromacs-2023.5/bin/gmx)
#   TOP_NDX - index file for trjconv/trjcat (default: index.ndx, in this dir)
#   NK      - which segment's .tpr to use as -s reference (default: 1)
#   SKIP    - frame skip (ps) for the VMD trajectory (default: 100)
#   CFG     - system_groups.cfg (for CONCAT_* settings)
#             (default: ../../0_prep/system_groups.cfg)
# =============================================================================

set -e

GROM="${GROM:-/gromacs/gromacs-2023.5/bin/gmx}"
NDX="${TOP_NDX:-index.ndx}"
NK="${NK:-1}"
SKIP="${SKIP:-100}"
CFG="${CFG:-../../0_prep/system_groups.cfg}"

log() { echo -e "[$(date '+%H:%M:%S')] $*"; }
die() { echo -e "[ERROR] $*" >&2; exit 1; }

[ -f "$NDX" ] || die "index file not found: $NDX"
[ -f "$CFG" ] || die "system_groups.cfg not found: $CFG"
command -v "$GROM" >/dev/null 2>&1 || die "gmx binary not found: $GROM"

source "$CFG"
: "${CONCAT_HAS_SOLUTE:=1}"
: "${CONCAT_SOLUTE_GROUP:=Protein}"
: "${CONCAT_CENTER_OUTPUT_GROUP:=Protein}"

NK_PAD=$(printf "%03d" "$NK")

# ── FIND SEGMENTS TO CONCATENATE ──────────────────────────────────────────────
# NOTE: 000.gro is reserved by the production script as a symlink to the
# starting structure (NPT endpoint) — it is created via "ln -s ... 000.gro"
# and is NEVER a real finished production segment, so it never has a
# matching 000.xtc. Real production output always starts at 001. Without
# excluding it here, find's glob "[0-9][0-9][0-9].gro" matches 000.gro too
# (find matches symlinks by name, not just real files), which makes this
# script mistake the starting structure for a finished segment and die
# looking for a 000.xtc that was never supposed to exist.
mapfile -t SEG_FILES < <(find . -maxdepth 1 -name '[0-9][0-9][0-9].gro' ! -name '000.gro' | sed 's#^\./##' | sort)
[ "${#SEG_FILES[@]}" -gt 0 ] || die "no segment .gro files found in $(pwd) (segments start at 001; 000.gro is the starting structure, not a segment)"
NB=${#SEG_FILES[@]}
log "Found $NB completed segment(s): ${SEG_FILES[*]}"

INPUT_XTCS=""
for f in "${SEG_FILES[@]}"; do
    seg="${f%.gro}"
    [ -f "${seg}.xtc" ] || die "expected trajectory ${seg}.xtc not found (segment ${seg} has a .gro but no .xtc)"
    INPUT_XTCS="${INPUT_XTCS}${seg}.xtc "
done

# ── FIND THE REFERENCE .tpr (handles gzip'd tprs — production compresses them) ─
TPR="${NK_PAD}.tpr"
CLEANUP_TPR=0
if [ ! -f "$TPR" ] && [ -f "${TPR}.gz" ]; then
    log "Decompressing ${TPR}.gz for use as -s reference..."
    gunzip -k "${TPR}.gz"
    CLEANUP_TPR=1
fi
[ -f "$TPR" ] || die "reference tpr not found: $TPR (or ${TPR}.gz)"

log "=========================================================="
log "Concatenating $NB segment(s) — reference tpr: $TPR"
log "=========================================================="

# ── A. CONCATENATE RAW SEGMENTS ────────────────────────────────────────────────
CONFIRM=$(for ((i=0; i<NB; i++)); do echo "c"; done)
echo "$CONFIRM" | "$GROM" trjcat -f $INPUT_XTCS -o traj_all_no_PBC.xtc -settime

# ── B. CENTER (IF SOLUTE) + APPLY PBC ─────────────────────────────────────────
if [ "$CONCAT_HAS_SOLUTE" -eq 0 ]; then
    "$GROM" trjconv -f traj_all_no_PBC.xtc -n "$NDX" -s "$TPR" -pbc mol \
        -o traj_all_center_PBC.xtc <<EOF
$CONCAT_CENTER_OUTPUT_GROUP
EOF
else
    "$GROM" trjconv -f traj_all_no_PBC.xtc -n "$NDX" -s "$TPR" -center -pbc mol \
        -o traj_all_center_PBC.xtc <<EOF
$CONCAT_SOLUTE_GROUP
$CONCAT_CENTER_OUTPUT_GROUP
EOF
fi

log "Removing intermediate traj_all_no_PBC.xtc..."
rm -f traj_all_no_PBC.xtc

# ── C. SKIP-REDUCED TRAJECTORY FOR VISUALIZATION ──────────────────────────────
"$GROM" trjconv -f traj_all_center_PBC.xtc -s "$TPR" -n "$NDX" -skip "$SKIP" \
    -o traj_all_center_PBC_VMD.xtc <<EOF
$CONCAT_CENTER_OUTPUT_GROUP
EOF

[ "$CLEANUP_TPR" -eq 1 ] && rm -f "$TPR"

echo
log "✅ Done."
log "   Analysis trajectory      : traj_all_center_PBC.xtc"
log "   Visualization trajectory : traj_all_center_PBC_VMD.xtc"