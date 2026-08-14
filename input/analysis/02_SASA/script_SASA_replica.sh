#!/bin/bash -e
# =============================================================================
# script_SASA_replica.sh
#
# Runs gmx sasa on the "Protein" group (GROMACS's built-in default group -
# no custom index needed) for the single replica this script is run from,
# then produces a SASA vs time plot + probability-density distribution via
# plot_sasa.py.
#
# Run from inside a replica's analysis subfolder, e.g. 4_prd/rep1/analysis/05_SASA/
#
# Env vars (all optional, fallback defaults match that layout):
#   GROM       - gmx binary                 (default: gmx)
#   PRD_DIR    - replica root directory     (default: ../../)
#   TPR        - .tpr file                  (default: auto-detected in PRD_DIR)
#   TRAJ       - trajectory file            (default: auto-detected in PRD_DIR)
#   REP_LABEL  - label for output filenames (default: basename of PRD_DIR)
#   OUT_DIR    - where outputs go           (default: current directory)
#   BIN_WIDTH  - override auto (Freedman-Diaconis) distribution bin width
# =============================================================================

set -e

GROM="${GROM:-gmx}"
NDX="${NDX:-../00_index/index_analysis.ndx}"
PRD_DIR="${PRD_DIR:-../..}"
OUT_DIR="${OUT_DIR:-.}"
GROUP="Protein"

log() { echo -e "[$(date '+%H:%M:%S')] $*"; }
die() { echo -e "[ERROR] $*" >&2; exit 1; }

mkdir -p "$OUT_DIR"
[ -d "$PRD_DIR" ] || die "replica directory not found: $PRD_DIR"
REP_LABEL="${REP_LABEL:-$(basename "$(cd "$PRD_DIR" && pwd)")}"

# ── AUTO-DETECT TPR ────────────────────────────────────────────────────────────
if [ -z "$TPR" ]; then
    TPR=$(find "$PRD_DIR" -maxdepth 1 -name "*.tpr" | head -n1)
fi
[ -n "$TPR" ] && [ -f "$TPR" ] || die "no .tpr found in $PRD_DIR — set TPR explicitly."

# ── LOCATE TRAJECTORY ──────────────────────────────────────────────────────────
# Same fixed name concat_trj_replica.sh always writes — no fallback to a raw
# segment .xtc, since that would silently run SASA on incomplete data with
# no warning (same convention as script_Rg_replica.sh).
if [ -z "$TRAJ" ]; then
    TRAJ="$PRD_DIR/traj_all_center_PBC.xtc"
fi
[ -f "$TRAJ" ] || die "concatenated trajectory not found: $TRAJ — has concat_trj_replica.sh run for this replica yet? (check for a .concat_done marker in $PRD_DIR, or set TRAJ explicitly)"

[ -f "$NDX" ] || die "index file not found: $NDX (run index_maker.sh first)"

log "==============================================="
log "SASA analysis — replica: $REP_LABEL"
log "==============================================="
log "TPR   : $TPR"
log "Traj  : $TRAJ"
log "Index : $NDX"
log "Group : $GROUP"
echo

outfile="$OUT_DIR/SASA_${GROUP}_${REP_LABEL}.xvg"
log "→ group '$GROUP' → $outfile"
echo "$GROUP" | "$GROM" sasa \
    -s "$TPR" \
    -n "$NDX" \
    -f "$TRAJ" \
    -o "$outfile"

log "✅ Finished SASA calculation for replica $REP_LABEL"

# ── PLOT + DISTRIBUTION ────────────────────────────────────────────────────────
PLOT_SCRIPT="$(dirname "$0")/plot_sasa.py"
if [ -f "$PLOT_SCRIPT" ]; then
    echo
    log "Generating SASA vs time + distribution plots..."
    BIN_ARGS=()
    [ -n "$BIN_WIDTH" ] && BIN_ARGS=(--bin-width "$BIN_WIDTH")
    python3 "$PLOT_SCRIPT" --outdir "$OUT_DIR" --label "$REP_LABEL" "${BIN_ARGS[@]}" "$outfile"
else
    log "WARNING: plot_sasa.py not found next to this script — skipping plot."
fi