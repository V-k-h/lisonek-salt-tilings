#!/bin/bash
# Sharded dim_m2_check.py on Trillium: each array task is one 192-core node running 192 shards, so
#   sbatch --array=0-(N-1) m2_dim4.sbatch P D
# runs the whole check for q = P^2, dim K = D as 192*N shards on N nodes.  Submit from $SCRATCH/ks-quadric
# ($HOME is read-only on the compute nodes).  Shards whose log already has a RESULT line are skipped, so
# resubmitting the same command after a timeout only reruns the unfinished ones.
# Merge afterwards:  python merge_m2.py logs/q${P}_d${D}/*_of_$((192*N)).log
#SBATCH --job-name=m2check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=24:00:00
#SBATCH --array=0
#SBATCH --output=logs/%x_%A_%a.out

P=${1:-5}; D=${2:-4}
TASK=${SLURM_ARRAY_TASK_ID:-0}; NSH=$((192 * ${SLURM_ARRAY_TASK_COUNT:-1}))
module load python/3.13
source $HOME/.virtualenvs/ks/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1   # one core per shard

mkdir -p logs/q${P}_d${D}
for c in $(seq 0 191); do
    i=$((TASK * 192 + c)); f=logs/q${P}_d${D}/shard_${i}_of_${NSH}.log
    grep -q '^RESULT' "$f" 2>/dev/null && continue
    python -u dim_m2_check.py $P $D $i $NSH > "$f" 2>&1 &
done
wait
echo "node $TASK: $(grep -l '^RESULT' logs/q${P}_d${D}/*_of_${NSH}.log | wc -l) of $NSH shards finished so far"
