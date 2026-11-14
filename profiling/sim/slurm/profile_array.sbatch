#!/bin/bash
# One shard of a sweep per array task. Submit with submit_array.sh, which also
# queues the merge; or directly:
#   BASE=generated/sim-profile/big sbatch --array=0-3 profiling/sim/slurm/profile_array.sbatch
# Datasets are split by series (every gate count of one qubit count / gate set),
# round-robin by qubit count, so auto-skip still works within each shard.
#SBATCH --job-name=sim-profile
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --exclusive
#SBATCH --signal=B:TERM@300
#SBATCH --output=slurm-%x-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source profiling/sim/slurm/env.sh

: "${BASE:?set BASE to the directory that holds the shard runs}"
SWEEP=${SWEEP:-profiling/sim/sweeps/scaling.toml}
BUILD_DIR=${BUILD_DIR:-build-slurm}
ID=${SLURM_ARRAY_TASK_ID:?run as an array job}
COUNT=${SLURM_ARRAY_TASK_COUNT:-$(( ${SLURM_ARRAY_TASK_MAX:?} + 1 ))}
RUN_DIR="$BASE/shard-$ID"
TOOLS=(--generator "$BUILD_DIR/src/generator/qcanopy_generate" --bench "$BUILD_DIR/profiling/sim/sim_bench")

if [[ -f "$RUN_DIR/run.json" ]]; then
  ARGS=(run --resume "$RUN_DIR" "${TOOLS[@]}")
else
  ARGS=(run --sweep "$SWEEP" --shard "$ID/$COUNT" --run-dir "$RUN_DIR" "${TOOLS[@]}" "$@")
fi

"$PYTHON" profiling/sim/sim_profile.py "${ARGS[@]}" &
PID=$!
trap 'kill -TERM "$PID" 2>/dev/null' TERM INT
set +e
wait "$PID"; STATUS=$?
if kill -0 "$PID" 2>/dev/null; then wait "$PID"; STATUS=$?; fi  # trap fired; wait for the runner
set -e
exit "$STATUS"
