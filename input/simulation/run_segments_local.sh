#!/bin/bash -e
# =============================================================================
# run_segments_local.sh — sequential (non-SLURM) segment runner
#
# Re-implements, for local execution, the SAME per-segment logic used by
# slurm_MD_eq_NPT_allatom.sh / slurm_MD_prd_allatom.sh:
#   - segment numbering is detected from existing NNN.gro files
#   - tinit is computed as nsteps * dt * (segment-1), read directly from
#     the .mdp file itself (same awk one-liner as the SLURM scripts)
#   - grompp + mdrun, continuing each segment from the previous one's .gro
#
# What is deliberately NOT reproduced here (SLURM-only concerns):
#   - sbatch submission / partitions / requeue
#   - staging to /tmp local scratch and copying back
#   - self-requeue (./script.sh calling itself) — replaced by a plain for-loop
#   - gzip of finished segment files
#
# This script does not alter the SLURM scripts or .dat files in any way and
# is only invoked when EXECUTION_MODE=local in run_setup_ab42.sh.
#
# Required environment variables:
#   GMX        gmx executable to use
#   NTHREADS   thread count for mdrun -nt
#   MDP        path to the (already nsteps/dt-patched) .mdp file to run
#   TOP        path to topology (topol.top)
#   NDX        path to index file (index.ndx)
#   GRO        starting structure to use for segment 000 (first run only)
#   SEGMENTS   number of segments (blocks) to run in total
# Optional:
#   MAXWARN    grompp -maxwarn value (default: 1)
#
# Run with the current working directory set to the target stage directory
# (e.g. 3_npt/, or 4_prd/repN/) — same convention as the SLURM scripts.
# =============================================================================

: "${GMX:?GMX not set}"
: "${NTHREADS:?NTHREADS not set}"
: "${MDP:?MDP not set}"
: "${TOP:?TOP not set}"
: "${NDX:?NDX not set}"
: "${GRO:?GRO not set}"
: "${SEGMENTS:?SEGMENTS not set}"
MAXWARN="${MAXWARN:-1}"

echo "[local] Starting sequential run: ${SEGMENTS} segment(s), mdp=${MDP}, gro0=${GRO}"

for ((i = 1; i <= SEGMENTS; i++)); do
    j=$(printf "%03d" "$i")
    k=$((i - 1))
    l=$(printf "%03d" "$k")

    # Segment 000 = the user-supplied starting structure (same convention
    # as the SLURM scripts: symlinked once, on the first segment only).
    if [ "$i" -eq 1 ] && [ ! -e 000.gro ]; then
        ln -s "$(readlink -f "$GRO")" 000.gro
    fi

    if [ -f "${j}.gro" ]; then
        echo "[local] Segment ${j} already present — skipping (resume behaviour)."
        continue
    fi

    echo "[local] ---------------------------------------------------------"
    echo "[local] Segment ${j} / $(printf "%03d" "$SEGMENTS")  (continuing from ${l}.gro)"

    # Same tinit calculation as the SLURM scripts: nsteps * dt * k, read
    # directly from the .mdp file so it always matches what will actually run.
    tinit=$(awk -v K="$k" '$1=="nsteps"{a=$3} $1=="dt"{b=$3} END{print a*b*K}' "$MDP")
    sed "s/tinit\s*=\s*0\.0/tinit = ${tinit}/g" "$MDP" > "${j}.mdp"

    "$GMX" grompp \
        -f  "${j}.mdp" \
        -po "${j}_out.mdp" \
        -c  "${l}" \
        -n  "$NDX" \
        -p  "$TOP" \
        -pp TMP_processed.top \
        -o  "${j}.tpr" \
        -maxwarn "$MAXWARN"

    "$GMX" mdrun \
        -nt   "$NTHREADS" \
        -pin  auto \
        -s    "${j}.tpr" \
        -x    "${j}.xtc" \
        -c    "${j}.gro" \
        -e    "${j}.edr" \
        -g    "${j}.log"

    echo "[local] Segment ${j} complete."
done

echo "[local] All ${SEGMENTS} segment(s) complete."
