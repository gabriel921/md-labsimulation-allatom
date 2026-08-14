#!/bin/bash -e
# =============================================================================
# index_maker.sh
#
# Builds an analysis index (chains / peptides / backbone-only groups) on top
# of an existing .ndx, driven entirely by system_groups.cfg - no atom ranges
# hardcoded here. Doesn't rely on deleting pre-existing default groups by
# number (that count varies per system/protein). Instead it probes
# make_ndx's own group listing to find out where new group numbering should
# start, and resolves the Backbone group by name.
#
# Run from inside a replica's analysis subfolder, e.g. 4_prd/rep1/analysis/00_index/
#
# Env vars (all optional, fallback defaults match that layout):
#   GRO     - structure file             (default: ../../000.gro)
#   NDX_IN  - existing index to extend   (default: ../../index.ndx)
#   NDX_OUT - output index file          (default: index_analysis.ndx)
#   GROM    - gmx binary                 (default: gmx)
#   CFG     - system_groups.cfg to use   (default: ../../../../0_prep/system_groups.cfg)
# =============================================================================

set -e

GRO="${GRO:-../../001.gro}"
NDX_IN="${NDX_IN:-../../index.ndx}"
NDX_OUT="${NDX_OUT:-index_analysis.ndx}"
GROM="${GROM:-gmx}"
CFG="${CFG:-../../../../0_prep/system_groups.cfg}"

log() { echo -e "[$(date '+%H:%M:%S')] $*"; }
die() { echo -e "[ERROR] $*" >&2; exit 1; }

[ -f "$GRO" ]    || die "structure file not found: $GRO"
[ -f "$NDX_IN" ] || die "input index not found: $NDX_IN"
[ -f "$CFG" ]    || die "system_groups.cfg not found: $CFG"
command -v "$GROM" >/dev/null 2>&1 || die "gmx binary not found: $GROM"

source "$CFG"

log "Building analysis index: $NDX_OUT"
log "  Structure : $GRO"
log "  Base ndx  : $NDX_IN"
log "  Config    : $CFG"

# ── PROBE PASS - list existing groups so we know where new numbering starts ──
PROBE_LOG=$(mktemp)
PROBE_NDX=$(mktemp --suffix=.ndx)
trap 'rm -f "$PROBE_LOG" "$PROBE_NDX"' EXIT

echo "q" | "$GROM" make_ndx -f "$GRO" -n "$NDX_IN" -o "$PROBE_NDX" > "$PROBE_LOG" 2>&1 || true

# Group listing lines look like:  "  4 Backbone             :   170 atoms"
mapfile -t EXISTING < <(grep -E '^ *[0-9]+ +[^:]+ *:' "$PROBE_LOG")
NEXT_IDX=${#EXISTING[@]}
[ "$NEXT_IDX" -gt 0 ] || die "could not parse any existing groups from make_ndx output — check $GROM/$GRO/$NDX_IN manually. Probe log: $PROBE_LOG"

resolve_group() {
    local name="$1"
    grep -E "^ *[0-9]+ +${name} *:" "$PROBE_LOG" | head -n1 | awk '{print $1}'
}

BACKBONE_IDX=""
if [ "${#BACKBONE_FOR[@]}" -gt 0 ]; then
    BACKBONE_IDX=$(resolve_group "$BACKBONE_GROUP_NAME")
    if [ -z "$BACKBONE_IDX" ]; then
        log "Available default groups:"
        printf '  %s\n' "${EXISTING[@]}" >&2
        die "BACKBONE_GROUP_NAME='$BACKBONE_GROUP_NAME' not found among default groups (see list above)."
    fi
fi

log "  Existing default groups : $NEXT_IDX"
[ -n "$BACKBONE_IDX" ] && log "  Backbone group resolved : $BACKBONE_IDX ($BACKBONE_GROUP_NAME)"
echo

# ── BUILD THE make_ndx COMMAND SEQUENCE ──────────────────────────────────────
declare -A GROUP_IDX
CMDS=()
ALL_MIN=""
ALL_MAX=""

add_range_group() {
    local name="$1" range="$2"
    local start="${range%-*}" end="${range#*-}"

    CMDS+=("a ${range}")
    CMDS+=("name ${NEXT_IDX} ${name}")
    GROUP_IDX["$name"]=$NEXT_IDX
    NEXT_IDX=$((NEXT_IDX + 1))

    if [ -z "$ALL_MIN" ] || [ "$start" -lt "$ALL_MIN" ]; then ALL_MIN=$start; fi
    if [ -z "$ALL_MAX" ] || [ "$end" -gt "$ALL_MAX" ]; then ALL_MAX=$end; fi
}

for def in "${CHAIN_DEFS[@]}"; do
    add_range_group "${def%%:*}" "${def#*:}"
done

for def in "${PEPTIDE_DEFS[@]}"; do
    add_range_group "${def%%:*}" "${def#*:}"
done

if [ -n "${ALL_CHAINS_NAME:-}" ] && [ "${#CHAIN_DEFS[@]}" -gt 0 ]; then
    CMDS+=("a ${ALL_MIN}-${ALL_MAX}")
    CMDS+=("name ${NEXT_IDX} ${ALL_CHAINS_NAME}")
    NEXT_IDX=$((NEXT_IDX + 1))
fi

for name in "${BACKBONE_FOR[@]}"; do
    idx="${GROUP_IDX[$name]:-}"
    if [ -z "$idx" ]; then
        log "WARNING: BACKBONE_FOR references unknown group '$name' — skipping."
        continue
    fi
    CMDS+=("${BACKBONE_IDX} & ${idx}")
    CMDS+=("name ${NEXT_IDX} ${name}_BB")
    NEXT_IDX=$((NEXT_IDX + 1))
done

CMDS+=("q")

log "Running make_ndx with ${#CMDS[@]} generated commands..."
printf '%s\n' "${CMDS[@]}" | "$GROM" make_ndx -f "$GRO" -n "$NDX_IN" -o "$NDX_OUT"

echo
log "✅ Index file written: $NDX_OUT"
