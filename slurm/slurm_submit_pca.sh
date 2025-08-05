#!/bin/bash
#SBATCH --job-name=pca_tets            # Job name
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=nikhil.pathak@praexisio.com.tw # Where to send mail
#SBATCH --ntasks=1                     # Number of tasks to run
#SBATCH --cpus-per-task=2              # threads
#SBATCH --mem=5gb                      # Job memory request
#SBATCH --output=pca_test_%j.log        # Standard output and error log
#SBATCH --account=YangLab              # The account 
#SBATCH --time=00:20:00                # Estimated Finish Time in hour:min:sec

# Record the start time
start_time=$(date +%s)

# Execute the MD job
# singularity exec -B /raid --nv /raid/YangLab/nik/dockingxopenmm.sif \
#     python simulations.py \
#     -o /raid/YangLab/nik/workbench/protein_MD/mdrun/1azm \
#     -i /raid/YangLab/nik/workbench/drdocking/drdock_valid/1azm/Protein \
#     -c 0.1 \
#     -g 0 \
#     -t 100 \
#     -f 100

singularity exec -B /raid \
    --nv /raid/YangLab/nik/dockingxopenmm.sif bash
    python do_PCA_and_clustering.py \
    --ps "protein and not (name H*)" \
    -n pca_test \
    --traj /raid/YangLab/nik/workdata/docking/viper/rps/MD/prod_10.dcd \
    --top /raid/YangLab/nik/workdata/docking/viper/rps/MD/model.prmtop \
    --os "protein and not (name H*)"

# Record the end time and running time
end_time=$(date +%s)
running_time=$((end_time - start_time))
hours=$(($running_time / 3600))
minutes=$((($running_time / 60) % 60))
seconds=$(($running_time % 60))

# Print the total running time
echo "Total running time: $hours hours $minutes minutes $seconds seconds"
