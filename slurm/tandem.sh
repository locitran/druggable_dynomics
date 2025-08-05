#!/bin/bash
#SBATCH --job-name=notebook              # Job name
#SBATCH --mail-type=BEGIN                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=quangloctrandinh1998vn@gmail.com  # Where to send mail
#SBATCH --ntasks=1                       # Run on a single CPU
#SBATCH --cpus-per-task=8               # Cores
#SBATCH --gres=gpu:1g.5gb:1             # Request GPU "generic resources"
#SBATCH --mem=25gb                       # Job memory request
#SBATCH --output=notebook.log            # Standard output and error log
#SBATCH --partition=COMPUTE1Q            # The partition that job submit to
#SBATCH --account=YangLab                # The account name

port=$(getAvailablePort)

# Port forward to the login node
/usr/bin/ssh -N -f -R $port:localhost:$port yang_loci@a100
conda activate tandem
cd /mnt/nas_1/YangLab/loci/druggable_dynomics
###################################################################################################
python -m notebook --no-browser --allow-root --port $port --NotebookApp.allow_remote_access=True