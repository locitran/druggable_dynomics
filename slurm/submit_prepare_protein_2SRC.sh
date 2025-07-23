#!/bin/bash
#SBATCH --job-name=2SRC_prepare_docking             # Job name
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                     # Number of tasks to run
#SBATCH --cpus-per-task=1              # threads
#SBATCH --mem=5gb                      # Job memory request
#SBATCH --output=2SRC%j.log        # Standard output and error log
#SBATCH --account=YangLab              # The account 
#SBATCH --time=00:20:00                # Estimated Finish Time
#SBATCH --partition=COMPUTE1Q          # Specify queue (max.time: 4hours)

basedir=/mnt/nas_1/YangLab/loci/druggable_dynomics/
singularity exec \
    --bind /raid \
    --bind $basedir:/md \
    --nv /raid/YangLab/nik/dockingxopenmm.sif \
    python /md/prepare_protein.py \
    -w /md/output/2SRC \
    -i /md/output/2SRC/model.pdb