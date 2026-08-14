#!/bin/bash -e
# =============================================================================
# eq_NVT_allatom.sh — NVT equilibration
# Can be run standalone or called from run_setup_ab42.sh via env variables.
# =============================================================================

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# When called from run_setup_ab42.sh these are set by the master script.
# When run standalone, set them here.
GMX=${GMX:-gmx_mpi}
NTHREADS=${NTHREADS:-8}
MDP=${MDP:-../mdp/eq_nvt_allatom.mdp}
GRO=${GRO:-../1_min/em.gro}
TOP=${TOP:-../0_prep/topol.top}
NDX=${NDX:-index.ndx}
export GMXLIB=${GMXLIB:-$HOME/GMXLIB}
# ── END CONFIGURATION ────────────────────────────────────────────────────────

echo "[eq_NVT] Starting NVT equilibration..."
echo "  GMX      : $GMX"
echo "  MDP      : $MDP"
echo "  Input    : $GRO"
echo "  Topology : $TOP"

# Generate default index file if not already present
if [ ! -f "$NDX" ]; then
    echo "[eq_NVT] Index file not found — generating default index..."
    echo "q" | $GMX make_ndx -f $GRO -o $NDX
fi

$GMX grompp \
    -f $MDP \
    -c $GRO \
    -r $GRO \
    -p $TOP \
    -n $NDX \
    -o nvt.tpr \
    -maxwarn 1

$GMX mdrun \
    -s      nvt.tpr \
    -deffnm nvt \
    -ntmpi  1 \
    -ntomp  $NTHREADS \
    -pin    auto

echo "[eq_NVT] Done. Output: nvt.gro"
