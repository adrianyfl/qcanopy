#!/bin/bash
# One sweep in one job. Resubmitting the same command resumes the same run.
#   sbatch profiling/sim/slurm/profile.sbatch [sim_profile.py run options...]
#   SWEEP=profiling/sim/sweeps/strategies.toml sbatch profiling/sim/slurm/profile.sbatch
#   sbatch --mem=200G profiling/sim/slurm/profile.sbatch --set 'generator.qubits=[26,28,30]'
# Env:  SWEEP, RUN_DIR (default generated/sim-profile/<sweep>-slurm), BUILD_DIR
#
# The simulator is single-threaded and memory-bound: 2 CPUs (runner + worker)
# are enough, memory decides the largest qubit count (see README), and
# --exclusive gives stable timings by keeping other jobs off the node's
# memory bus (drop it if your queue is long).
#SBATCH --job-name=sim-profile
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --exclusive
#SBATCH --signal=B:TERM@300
#SBATCH --output=slurm-%x-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source profiling/sim/slurm/env.sh

SWEEP=${SWEEP:-profiling/sim/sweeps/scaling.toml}
RUN_DIR=${RUN_DIR:-generated/sim-profile/$(basename "$SWEEP" .toml)-slurm}
BUILD_DIR=${BUILD_DIR:-build-slurm}
TOOLS=(--generator "$BUILD_DIR/src/generator/qcanopy_generate" --bench "$BUILD_DIR/profiling/sim/sim_bench")

if [[ -f "$RUN_DIR/run.json" ]]; then
  echo "resuming $RUN_DIR"
  ARGS=(run --resume "$RUN_DIR" "${TOOLS[@]}")
else
  ARGS=(run --sweep "$SWEEP" --run-dir "$RUN_DIR" "${TOOLS[@]}" "$@")
fi

# Slurm signals the batch shell 5 minutes before the time limit; pass it on so
# the runner stops its workers and leaves the run resumable.
"$PYTHON" profiling/sim/sim_profile.py "${ARGS[@]}" &
PID=$!
trap 'kill -TERM "$PID" 2>/dev/null' TERM INT
set +e
wait "$PID"; STATUS=$?
if kill -0 "$PID" 2>/dev/null; then wait "$PID"; STATUS=$?; fi  # trap fired; wait for the runner
set -e
if [[ $STATUS -ne 0 ]]; then
  echo "stopped (exit ${STATUS}); resubmit the same command to resume $RUN_DIR" >&2
  exit "$STATUS"
fi
