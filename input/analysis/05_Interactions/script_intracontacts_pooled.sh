#!/bin/bash -e
# =============================================================================
# script_intracontacts_pooled.sh
#
# Thin launcher for the pooled intra-chain contacts analysis. Auto-discovers
# every replica's production trajectory under 4_prd/, points at the shared
# nowater.tpr / index_analysis.ndx, and hands everything to
# intra_contacts_pooled.py via environment variables.
#
# Expects to be run from its working directory (4_prd/Analysis/<subfolder>/),
# as dispatched by run_setup_analysis.sh — same convention as the per-replica
# analysis scripts and script_hdbscan_pooled.sh.
#
# Env vars accepted (all optional, sensible defaults shown):
#   PRD_DIR       - path to 4_prd/ (default: computed as ../.. from cwd)
#   TRAJ_NAME     - per-replica trajectory filename (default: traj_all_center_PBC.xtc)
#   NDXFILE       - index file (default: $PRD_DIR/rep1/analysis/00_index/index_analysis.ndx)
#   TOPOLOGY      - MDAnalysis topology (default: $PRD_DIR/Analysis/nowater.tpr)
#   CHAIN_GROUPS  - comma-separated NDX group names to analyze independently
#                   (default: Protein). Set this if you have multiple chains,
#                   e.g. CHAIN_GROUPS="ChainA,ChainB,ChainC,ChainD"
#   OUTDIR        - output subfolder for this script's results (default: contacts_intra)
#   PYTHON_BIN    - python3 interpreter to use (default: python3)
#
# All other intra_contacts_pooled.py options (CUTOFF, RESIDUE_EXCLUDE,
# ANALYSIS_INTERVAL_NS, MA_WINDOW, etc.) can also be exported before calling
# this script — they pass through untouched.
# =============================================================================

set -e

log() { echo -e "[$(date '+%H:%M:%S')] [intracontacts_pooled] $*"; }
die() { echo -e "[ERROR] [intracontacts_pooled] $*" >&2; exit 1; }

PRD_DIR="${PRD_DIR:-$(cd ../.. && pwd)}"
TRAJ_NAME="${TRAJ_NAME:-traj_all_center_PBC.xtc}"
NDXFILE="${NDXFILE:-$PRD_DIR/rep1/analysis/00_index/index_analysis.ndx}"
TOPOLOGY="${TOPOLOGY:-$PRD_DIR/Analysis/nowater.tpr}"
CHAIN_GROUPS="${CHAIN_GROUPS:-Protein}"
OUTDIR="${OUTDIR:-contacts_intra}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

[ -d "$PRD_DIR" ] || die "PRD_DIR not found: $PRD_DIR"
[ -f "$NDXFILE" ] || die "index file not found: $NDXFILE"
[ -f "$TOPOLOGY" ] || die "nowater topology not found: $TOPOLOGY (did the nowater prep step run?)"

# ── DISCOVER REPLICAS ─────────────────────────────────────────────────────
mapfile -t REPLICAS < <(find "$PRD_DIR" -maxdepth 1 -mindepth 1 -type d -name 'rep*' -exec basename {} \; | sort)
[ "${#REPLICAS[@]}" -gt 0 ] || die "no replicate folders found under $PRD_DIR"
log "Replicas found: ${REPLICAS[*]}"

# ── BUILD TRAJECTORY LIST (hard die on any missing trajectory — no silent
#    fallback to whatever .xtc happens to be alphabetically first) ─────────
TRAJFILES=""
REPNAMES=""
for rep in "${REPLICAS[@]}"; do
    tfile="$PRD_DIR/$rep/$TRAJ_NAME"
    [ -f "$tfile" ] || die "expected trajectory not found: $tfile (TRAJ_NAME=$TRAJ_NAME)"
    if [ -z "$TRAJFILES" ]; then
        TRAJFILES="$tfile"
        REPNAMES="$rep"
    else
        TRAJFILES="$TRAJFILES,$tfile"
        REPNAMES="$REPNAMES,$rep"
    fi
done

log "Topology     : $TOPOLOGY"
log "Index        : $NDXFILE"
log "Chain groups : $CHAIN_GROUPS"
log "Output to    : $OUTDIR"
log "Trajectories:"
IFS=',' read -ra _TF <<< "$TRAJFILES"
for t in "${_TF[@]}"; do log "  - $t"; done

export TOPOLOGY NDXFILE TRAJFILES REPNAMES CHAIN_GROUPS OUTDIR

"$PYTHON_BIN" intra_contacts_pooled.py