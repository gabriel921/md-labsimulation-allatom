#!/bin/bash -e
# =============================================================================
# min_allatom.sh — Energy minimisation
# Can be run standalone or called from run_setup_ab42.sh via env variables.
# =============================================================================

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# When called from run_setup_ab42.sh these are set by the master script.
# When run standalone, set them here.
GMX=${GMX:-gmx_mpi}
NTHREADS=${NTHREADS:-8}
MDP=${MDP:-../mdp/min_allatom.mdp}
GRO=${GRO:-../0_prep/ab42_solv_ions.gro}
TOP=${TOP:-../0_prep/topol.top}
export GMXLIB=${GMXLIB:-$HOME/GMXLIB}
# ── END CONFIGURATION ────────────────────────────────────────────────────────

echo "[min_allatom] Starting energy minimisation..."
echo "  GMX      : $GMX"
echo "  MDP      : $MDP"
echo "  Input    : $GRO"
echo "  Topology : $TOP"

$GMX grompp \
    -f $MDP \
    -c $GRO \
    -p $TOP \
    -o em.tpr \
    -maxwarn 1

$GMX mdrun \
    -s      em.tpr \
    -deffnm em \
    -ntmpi  1 \
    -ntomp  $NTHREADS \
    -pin    auto

echo "[min_allatom] Done. Output: em.gro"
