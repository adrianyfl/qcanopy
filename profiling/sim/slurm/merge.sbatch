#!/bin/bash
# Merge the shard runs under $BASE into $BASE/merged (summary + report.html).
#SBATCH --job-name=sim-profile-merge
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=slurm-%x-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source profiling/sim/slurm/env.sh
: "${BASE:?set BASE}"
rm -rf "$BASE/merged"
"$PYTHON" profiling/sim/sim_profile.py merge "$BASE/merged" "$BASE"/shard-*
