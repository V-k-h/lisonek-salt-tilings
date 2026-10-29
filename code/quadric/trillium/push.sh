#!/bin/bash
# Run on the laptop: copy the code and the Trillium scripts to $SCRATCH/ks-quadric on Trillium.
# Needs the "trillium" host entry in ~/.ssh/config (with ControlMaster, the MFA prompt comes only once).
set -e
cd "$(dirname "$0")/.."
R=$(ssh trillium 'echo $SCRATCH')/ks-quadric
ssh trillium "mkdir -p $R/logs"
rsync -av dim_m2_check.py merge_m2.py cross_tiling.py quadric_ks.py quadric_gh.py fq.py \
      trillium/m2_dim4.sbatch trillium/setup_trillium.sh trillium:$R/
echo "copied to trillium:$R"
