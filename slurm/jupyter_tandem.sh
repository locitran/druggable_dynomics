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

# Get an available port
port=$(getAvailablePort)

# Port forward to the login node
/usr/bin/ssh -N -f -R $port:localhost:$port yang_loci@a100

# notebook_link -i notebook.log

conda activate tandem
cd /mnt/nas_1/YangLab/loci/druggable_dynomics
# python -m notebook --no-browser --allow-root --port $port --NotebookApp.allow_remote_access=True
jupyter notebook --no-browser --allow-root --port $port --NotebookApp.allow_remote_access=True



jupyter notebook --ip=0.0.0.0 --no-browser

echo $port


http://localhost:8888/tree?token=5f97fcf0dbb103b08769a677bc418e77dc0f22eb58a39e2f
http://127.0.0.1:8888/tree?token=5f97fcf0dbb103b08769a677bc418e77dc0f22eb58a39e2f

http://140.114.97.192:15457/tree?token=261e1ea358f4ba5f50d51cbb5c5365ed84845e274b120bd6

http://140.114.97.192:15457/tree?token=261e1ea358f4ba5f50d51cbb5c5365ed84845e274b120bd6


      http://a100:8888/tree?token=385e888b0b96bb97bc9b8f743ec9ac6a8da95cb13645ae73
        http://0.0.0.0:8888/tree?token=385e888b0b96bb97bc9b8f743ec9ac6a8da95cb13645ae73