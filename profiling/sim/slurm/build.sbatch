#!/bin/bash
# Build qcanopy_generate and sim_bench on a compute node, in Release.
#   sbatch profiling/sim/slurm/build.sbatch                 (from the repo root)
# Env:  BUILD_DIR (default build-slurm), CMAKE_ARGS (extra configure flags)
#SBATCH --job-name=qcanopy-build
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%x-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source profiling/sim/slurm/env.sh

BUILD_DIR=${BUILD_DIR:-build-slurm}
# shellcheck disable=SC2086
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release ${CMAKE_ARGS:-}
cmake --build "$BUILD_DIR" -j "${SLURM_CPUS_PER_TASK:-4}"
ls -l "$BUILD_DIR"/src/generator/qcanopy_generate "$BUILD_DIR"/profiling/sim/sim_bench
