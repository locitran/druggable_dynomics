import os
import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDWriter

# Define the path where the .dcd files are located
path = "/raid/YangLab/nik/workdata/mdrun/mtmc1_conf/mtmc1_c2b/prx_world_04797_14/"
output_file = os.path.join(path, 'merged_1_10.dcd')

# Initialize the universe with the topology file
prmfile = os.path.join(path, 'model.prmtop')
u = mda.Universe(prmfile)

# List to store trajectory filenames
trajectory_files = [os.path.join(path, f'prod_{i}.dcd') for i in range(1, 11)]

# Verify all trajectory files exist
for traj in trajectory_files:
    if not os.path.exists(traj):
        raise FileNotFoundError(f"File {traj} does not exist")

# Create a DCD writer
with DCDWriter(output_file, u.atoms.n_atoms) as writer:
    for traj in trajectory_files:
        # Load each trajectory and write frames to the output file
        u.load_new(traj)
        for ts in u.trajectory:
            writer.write(u)

print(f"Successfully merged trajectories into {output_file}")
