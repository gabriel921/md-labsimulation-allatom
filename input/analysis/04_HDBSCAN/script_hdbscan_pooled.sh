#!/bin/bash -e
# =============================================================================
# script_hdbscan_pooled.sh
#
# Thin launcher for the pooled HDBSCAN analysis. Auto-discovers every
# replica's production trajectory under 4_prd/, points at the shared
# nowater.gro / index_analysis.ndx, and hands everything to
# hdbscan_pooled.py via environment variables.
#
# Expects to be run from its working directory (4_prd/Analysis/<subfolder>/),
# as dispatched by run_setup_analysis.sh — same convention as the per-replica
# analysis scripts.
#
# Env vars accepted (all optional, sensible defaults shown):
#   PRD_DIR      - path to 4_prd/ (default: computed as ../.. from cwd)
#   TRAJ_NAME    - per-replica trajectory filename (default: traj_all_center_PBC.xtc)
#   NDXFILE      - index file (default: $PRD_DIR/rep1/analysis/00_index/index_analysis.ndx)
#   TOPFILE      - nowater topology (default: $PRD_DIR/Analysis/nowater.gro)
#   OUTDIR       - output subfolder for this script's results (default: reps_monomer)
#   PYTHON_BIN   - python3 interpreter to use (default: python3)
#
# All other hdbscan_pooled.py options (ATOM_NAME_FILTER, PCA_COMPONENTS,
# STRIDE, FRAME_START/END, etc.) can also be exported before calling this
# script — they pass through untouched.
# =============================================================================

set -e

log() { echo -e "[$(date '+%H:%M:%S')] [hdbscan_pooled] $*"; }
die() { echo -e "[ERROR] [hdbscan_pooled] $*" >&2; exit 1; }

PRD_DIR="${PRD_DIR:-$(cd ../.. && pwd)}"
TRAJ_NAME="${TRAJ_NAME:-traj_all_center_PBC.xtc}" #depending on the amount of frames in the trajectory you might want to change this to the VMD traj 
NDXFILE="${NDXFILE:-$PRD_DIR/rep1/analysis/00_index/index_analysis.ndx}"
TOPFILE="${TOPFILE:-$PRD_DIR/Analysis/nowater.gro}"
OUTDIR="${OUTDIR:-reps_monomer}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

[ -d "$PRD_DIR" ] || die "PRD_DIR not found: $PRD_DIR"
[ -f "$NDXFILE" ] || die "index file not found: $NDXFILE"
[ -f "$TOPFILE" ] || die "nowater topology not found: $TOPFILE (did the nowater prep step run?)"

# ── DISCOVER REPLICAS ─────────────────────────────────────────────────────
mapfile -t REPLICAS < <(find "$PRD_DIR" -maxdepth 1 -mindepth 1 -type d -name 'rep*' -exec basename {} \; | sort)
[ "${#REPLICAS[@]}" -gt 0 ] || die "no replicate folders found under $PRD_DIR"
log "Replicas found: ${REPLICAS[*]}"

# ── BUILD TRAJECTORY LIST (hard die on any missing trajectory — no silent
#    fallback to whatever .xtc happens to be alphabetically first) ─────────
TRAJFILES=""
for rep in "${REPLICAS[@]}"; do
    tfile="$PRD_DIR/$rep/$TRAJ_NAME"
    [ -f "$tfile" ] || die "expected trajectory not found: $tfile (TRAJ_NAME=$TRAJ_NAME)"
    if [ -z "$TRAJFILES" ]; then
        TRAJFILES="$tfile"
    else
        TRAJFILES="$TRAJFILES,$tfile"
    fi
done

log "Topology  : $TOPFILE"
log "Index     : $NDXFILE"
log "Output to : $OUTDIR"
log "Trajectories:"
IFS=',' read -ra _TF <<< "$TRAJFILES"
for t in "${_TF[@]}"; do log "  - $t"; done

export TOPFILE NDXFILE TRAJFILES OUTDIR

"$PYTHON_BIN" hdbscan_pooled.py