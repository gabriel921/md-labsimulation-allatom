#!/bin/bash -e
# =============================================================================
# prep_nowater.sh
#
# One-time (per system) prep step for POOLED analyses (HDBSCAN, intra-chain
# contacts, etc). Produces a water-stripped topology + structure that those
# analyses load once and reuse across all replicas:
#
#   $PRD_DIR/Analysis/nowater.tpr
#   $PRD_DIR/Analysis/nowater.gro
#
# Source files are pulled from rep1 only — topology is identical across
# replicas, so there's no need to regenerate this per replica.
#
# Idempotent: skipped entirely if both outputs already exist, unless
# FORCE_REGEN=1 is set.
#
# Expects to be run with these already set/exported by the caller
# (run_setup_analysis.sh):
#   PRD_DIR        - path to 4_prd/
#   GROM           - gmx binary to use (pinned to 2021.2 by convention,
#                    same as DSSP, unless overridden by the caller)
# =============================================================================

set -e

PRD_DIR="${PRD_DIR:?PRD_DIR must be set by the caller}"
GROM="${GROM:?GROM must be set by the caller}"
FORCE_REGEN="${FORCE_REGEN:-0}"

log() { echo -e "[$(date '+%H:%M:%S')] [prep_nowater] $*"; }
die() { echo -e "[ERROR] [prep_nowater] $*" >&2; exit 1; }

SRC_REP_DIR="$PRD_DIR/rep1"
SRC_TPR="$SRC_REP_DIR/001.tpr"
SRC_GRO="$SRC_REP_DIR/001.gro"
SRC_NDX="$SRC_REP_DIR/index.ndx"

OUT_DIR="$PRD_DIR/Analysis"
OUT_TPR="$OUT_DIR/nowater.tpr"
OUT_GRO="$OUT_DIR/nowater.gro"

mkdir -p "$OUT_DIR"

if [ "$FORCE_REGEN" -ne 1 ] && [ -f "$OUT_TPR" ] && [ -f "$OUT_GRO" ]; then
    log "nowater.tpr and nowater.gro already present in $OUT_DIR — skipping (FORCE_REGEN=0)."
    exit 0
fi

[ -f "$SRC_TPR" ] || die "source tpr not found: $SRC_TPR"
[ -f "$SRC_GRO" ] || die "source gro not found: $SRC_GRO"
[ -f "$SRC_NDX" ] || die "index not found: $SRC_NDX (expected per-replica index build to have run already)"

log "Generating nowater topology/structure from $SRC_REP_DIR (group: Protein)"

log "  convert-tpr..."
echo "Protein" | "$GROM" convert-tpr -s "$SRC_TPR" -n "$SRC_NDX" -o "$OUT_TPR" \
    || die "gmx convert-tpr failed"

log "  trjconv..."
echo "Protein" | "$GROM" trjconv -f "$SRC_GRO" -s "$SRC_TPR" -n "$SRC_NDX" -o "$OUT_GRO" \
    || die "gmx trjconv failed"

[ -f "$OUT_TPR" ] || die "convert-tpr reported success but $OUT_TPR is missing"
[ -f "$OUT_GRO" ] || die "trjconv reported success but $OUT_GRO is missing"

log "Done. → $OUT_TPR"
log "Done. → $OUT_GRO"