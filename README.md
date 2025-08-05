This repository contains the source code for the Druggable DynOmics project.

Github repository: https://github.com/locitran/druggable_dynomics.git

# 1. Run image/container interactive mode

```bash
basedir=/mnt/nas_1/YangLab/loci/druggable_dynomics/
singularity exec \
    --home /tmp \
    --bind /raid \
    --bind $basedir:/md \
    --nv /raid/images/dockingxopenmm \
    bash -c "cd /md && exec bash"
```

# 2. Writable mode (to install new packages)

**For example**: Install notebook

Run sandbox with `--writable` mode:
```bash
singularity exec \
    --home /tmp \
    --writable \
    --nv /raid/images/dockingxopenmm bash
```

Get out of the sandbox to run as interactive mode:
```bash
basedir=/mnt/nas_1/YangLab/loci/druggable_dynomics/
singularity exec \
    --home /tmp \
    --bind /raid \
    --bind $basedir:/md \
    --nv /raid/images/dockingxopenmm bash

# Test notebook
python -m notebook
```

# 3. Run jupyter notebook using `dockingxopenmm` sandbox

```bash
#!/bin/bash
#SBATCH --ntasks=1                       # Run on a single CPU
#SBATCH --cpus-per-task=8               # Cores
#SBATCH --gres=gpu:1g.5gb:1             # Request GPU "generic resources"
#SBATCH --mem=25gb                       # Job memory request
#SBATCH --output=notebook.log            # Standard output and error log

# Get an available port
port=$(getAvailablePort)

# Port forward to the login node
/usr/bin/ssh -N -f -R $port:localhost:$port yang_loci@a100

# Let's go
basedir=/mnt/nas_1/YangLab/loci/druggable_dynomics/
singularity exec \
    --home /tmp \
    --bind /raid \
    --bind $basedir:/md \
    --nv /raid/images/dockingxopenmm \
    bash -c "cd /md && python -m notebook --no-browser --allow-root --port $port --NotebookApp.allow_remote_access=True"

```

**Call the notebook**
```bash
notebook_link -i notebook.log
# To access the notebook,  copy and paste one of these URLs in a browser:
#          http://140.114.97.192:19419/?token=1a9022e70f8807b21dfa9dec88280f44ca06a43e88108730
```


# Useful `Git` command

```bash
git add .
git commit -m "msg"
git push origin main
```
