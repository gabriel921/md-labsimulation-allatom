#!/bin/bash -e
# =============================================================================
# script_Rg_replica.sh
#
# Runs gmx gyrate for every "*_BB" group found in the analysis index, for
# the single replica this script is run from, then produces an Rg vs time
# plot. No hardcoded group list - driven entirely by whatever groups
# index_maker.sh created from system_groups.cfg.
#
# Run from inside a replica's analysis subfolder, e.g. 4_prd/rep1/analysis/04_Rg/
#
# Env vars (all optional, fallback defaults match that layout):
#   GROM      - gmx binary                 (default: gmx)
#   NDX       - analysis index file        (default: ../00_index/index_analysis.ndx)
#   PRD_DIR   - replica root directory     (default: ../../)
#   TPR       - .tpr file                  (default: auto-detected in PRD_DIR)
#   TRAJ      - trajectory file            (default: auto-detected in PRD_DIR)
#   REP_LABEL - label for output filenames (default: basename of PRD_DIR)
#   OUT_DIR   - where outputs go           (default: current directory)
# =============================================================================

set -e

GROM="${GROM:-gmx}"
NDX="${NDX:-../00_index/index_analysis.ndx}"
PRD_DIR="${PRD_DIR:-../..}"
OUT_DIR="${OUT_DIR:-.}"

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
# concat_trj_replica.sh always writes this fixed name once all production
# segments are concatenated — no fallback to a raw segment .xtc, since that
# would silently run Rg on incomplete data with no warning.
if [ -z "$TRAJ" ]; then
    TRAJ="$PRD_DIR/traj_all_center_PBC.xtc"
fi
[ -f "$TRAJ" ] || die "concatenated trajectory not found: $TRAJ — has concat_trj_replica.sh run for this replica yet? (check for a .concat_done marker in $PRD_DIR, or set TRAJ explicitly)"

[ -f "$NDX" ] || die "index file not found: $NDX (run index_maker.sh first)"

log "==============================================="
log "Rg analysis — replica: $REP_LABEL"
log "==============================================="
log "TPR   : $TPR"
log "Traj  : $TRAJ"
log "Index : $NDX"
echo

# ── DISCOVER GROUPS ────────────────────────────────────────────────────────────
# NOTE: deliberately NOT named "GROUPS" — that's a reserved bash special
# variable (list of the current user's group memberships) and silently
# breaks array assignment if reused.
mapfile -t BB_GROUPS < <(grep -oE '^\[ *[^]]*_BB *\]' "$NDX" | tr -d '[]' | tr -d ' ')
[ "${#BB_GROUPS[@]}" -gt 0 ] || die "no *_BB groups found in $NDX — nothing to analyse."
log "Groups found: ${BB_GROUPS[*]}"
echo

XVG_FILES=()
for group in "${BB_GROUPS[@]}"; do
    outfile="$OUT_DIR/Rg_${group}_${REP_LABEL}.xvg"
    log "→ group '$group' → $outfile"
    echo "$group" | "$GROM" gyrate \
        -s "$TPR" \
        -n "$NDX" \
        -f "$TRAJ" \
        -o "$outfile"
    XVG_FILES+=("$outfile")
done

log "✅ Finished Rg calculation for replica $REP_LABEL"

# ── PLOT ──────────────────────────────────────────────────────────────────────
PLOT_SCRIPT="$(dirname "$0")/plot_rg.py"
if [ -f "$PLOT_SCRIPT" ]; then
    echo
    log "Generating Rg vs time plot..."
    python3 "$PLOT_SCRIPT" --outdir "$OUT_DIR" --label "$REP_LABEL" "${XVG_FILES[@]}"
else
    log "WARNING: plot_rg.py not found next to this script — skipping plot."
fi
