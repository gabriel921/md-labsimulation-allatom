#!/bin/bash -e
# =============================================================================
# script_dssp_replica.sh
#
# Runs gmx do_dssp (GROMACS 2021.2 syntax) on the "Protein" group for the
# single replica this script is run from, then hands the resulting .xpm
# matrix to plot_dssp.py for a per-residue secondary-structure distribution
# (and a bonus residue-vs-time heatmap).
#
# NOTE on GROMACS version: do_dssp's argument syntax changed between the
# GROMACS 2021.x line and later gmx dssp rewrites, so this deliberately
# calls a pinned 2021.2 gmx binary rather than "gmx" off PATH — override
# GROM if your 2021.2 install lives elsewhere.
#
# NOTE on the DSSP binary: do_dssp itself does not implement the DSSP
# algorithm — it shells out to an external `dssp`/`mkdssp` executable. Set
# DSSP_BIN below (or as an env var) if it isn't already on PATH/pointed to
# by a DSSP env var on the target machine.
#
# Run from inside a replica's analysis subfolder, e.g. 4_prd/rep1/analysis/06_DSSP/
#
# Env vars (all optional, fallback defaults match that layout):
#   GROM       - gmx binary (2021.2)        (default: /gromacs/gromacs-2021.2/bin/gmx)
#   DSSP_BIN   - external dssp/mkdssp binary (default: unset — assumes DSSP env
#                var or PATH is already configured on the target machine)
#   PRD_DIR    - replica root directory     (default: ../../)
#   GRO        - .GRO file                  (default: auto-detected in PRD_DIR)
#   TRAJ       - trajectory file            (default: auto-detected in PRD_DIR)
#   REP_LABEL  - label for output filenames (default: basename of PRD_DIR)
#   OUT_DIR    - where outputs go           (default: current directory)
# =============================================================================

set -e

export GMXLIB=$HOME/Desktop/GMXLIB

GROM="${GROM:-/gromacs/gromacs-2021.2/bin/gmx}"
NDX="${NDX:-../00_index/index_analysis.ndx}"
PRD_DIR="${PRD_DIR:-../..}"
OUT_DIR="${OUT_DIR:-.}"
GROUP="Protein"

log() { echo -e "[$(date '+%H:%M:%S')] $*"; }
die() { echo -e "[ERROR] $*" >&2; exit 1; }

[ -n "$DSSP_BIN" ] && export DSSP="$DSSP_BIN"

mkdir -p "$OUT_DIR"
[ -d "$PRD_DIR" ] || die "replica directory not found: $PRD_DIR"
REP_LABEL="${REP_LABEL:-$(basename "$(cd "$PRD_DIR" && pwd)")}"

# ── AUTO-DETECT GRO ────────────────────────────────────────────────────────────
if [ -z "$GRO" ]; then
    GRO=$(find "$PRD_DIR" -maxdepth 1 -name "*.gro" | head -n1)
fi
[ -n "$GRO" ] && [ -f "$GRO" ] || die "no .gro found in $PRD_DIR — set GRO explicitly."

# ── LOCATE TRAJECTORY ──────────────────────────────────────────────────────────
# Same fixed name concat_trj_replica.sh always writes — no fallback to a raw
# segment .xtc, since that would silently run DSSP on incomplete data.
if [ -z "$TRAJ" ]; then
    TRAJ="$PRD_DIR/traj_all_center_PBC.xtc"
fi
[ -f "$TRAJ" ] || die "concatenated trajectory not found: $TRAJ — has concat_trj_replica.sh run for this replica yet? (check for a .concat_done marker in $PRD_DIR, or set TRAJ explicitly)"

log "==============================================="
log "DSSP analysis — replica: $REP_LABEL"
log "==============================================="
log "GROM  : $GROM"
log "GRO   : $GRO"
log "Traj  : $TRAJ"
log "Group : $GROUP (GROMACS built-in default — no custom index needed)"
[ -n "$DSSP" ] && log "DSSP  : $DSSP"
echo

xpm_out="$OUT_DIR/DSSP_${GROUP}_${REP_LABEL}.xpm"
xvg_out="$OUT_DIR/DSSP_${GROUP}_${REP_LABEL}.xvg"
dat_out="$OUT_DIR/DSSP_${GROUP}_${REP_LABEL}.dat"

log "→ group '$GROUP'"
log "    xpm : $xpm_out"
log "    xvg : $xvg_out"
log "    dat : $dat_out"

echo "$GROUP" | "$GROM" do_dssp \
    -s "$GRO" \
    -n "$NDX"\
    -f "$TRAJ" \
    -o "$xpm_out" \
    -sc "$xvg_out" \
    -ssdump "$dat_out"

log "✅ Finished DSSP calculation for replica $REP_LABEL"

# ── PER-RESIDUE DISTRIBUTION + HEATMAP ─────────────────────────────────────────
PLOT_SCRIPT="$(dirname "$0")/plot_dssp.py"
if [ -f "$PLOT_SCRIPT" ]; then
    echo
    log "Generating per-residue SS distribution"
    python3 "$PLOT_SCRIPT"
else
    log "WARNING: plot_dssp.py not found next to this script — skipping plot."
fi
