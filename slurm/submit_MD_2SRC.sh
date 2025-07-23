#!/bin/bash
#SBATCH --job-name=2SRC_100ns_md         # Job name
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                     # Number of tasks to run
#SBATCH --cpus-per-task=4              # threads
#SBATCH --mem=5gb                      # Job memory request
#SBATCH --gres=gpu:7g.40gb:1           # gputhreads
#SBATCH --output=2SRC_100ns_md_%j.log        # Standard output and error log
#SBATCH --account=YangLab              # The account 
#SBATCH --time=24:00:00                # Estimated Finish Time in hour:min:sec

# Record the start time
start_time=$(date +%s)
basedir=/mnt/nas_1/YangLab/loci/druggable_dynomics
# Execute the MD job
singularity exec \
    -B /raid \
    --bind $basedir:/md \
    --nv /raid/YangLab/nik/dockingxopenmm.sif \
    python /md/simulations.py \
    -o /md/output/2SRC \
    -i /md/output/2SRC/Protein \
    -g 0 \
    -c 0.1 \
    -t 100 \
    -f 100

# Record the end time and running time
end_time=$(date +%s)
running_time=$((end_time - start_time))
hours=$(($running_time / 3600))
minutes=$((($running_time / 60) % 60))
seconds=$(($running_time % 60))

# Print the total running time
echo "Total running time: $hours hours $minutes minutes $seconds seconds"