# File: keep_running_production.py
import os
import sys
from parmed import load_file
from openmm import *
from openmm.app import *
from openmm.unit import *
from parmed.openmm import RestartReporter
from OpenMMUtilities import *

# Parse command-line arguments
if len(sys.argv) < 5:
    print("Usage: python keep_running_production.py <path> <prmfile> <ncrst_file> <sim_indices>")
    print("Example: python keep_running_production.py /path/to/files model.prmtop prod_1.ncrst 2,3,4,5")
    sys.exit(1)

# Input parameters
path = sys.argv[1]
print(path)
prmfile = sys.argv[2]
ncrst_file = sys.argv[3]
sims = list(map(int, sys.argv[4].split(',')))

# System setup
prmtop = AmberPrmtopFile(prmfile)
#prmtop = AmberPrmtopFile(os.path.join(path, prmfile))
print(prmtop)
system = prmtop.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=10 * angstroms,
    constraints=HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005,
)

# Simulation parameters
dt = 1.0 * femtoseconds
fric_coeff = 1.0 / picosecond
temperature = 310.0 * kelvin
integrator = LangevinIntegrator(temperature, fric_coeff, dt)
prod_time = 100 * nanoseconds #time of each run
advanced_steps = int(prod_time / dt)
frames = 2 * advanced_steps

# Platform settings
platform = Platform.getPlatformByName('CUDA')
platformProp = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}

# Initialize simulation
simulation = Simulation(prmtop.topology, system, integrator, platform, platformProp)
simulation.context.getIntegrator().setTemperature(temperature)

# Load initial restart file
rst = load_file(ncrst_file)
#rst = load_file(os.path.join(path, ncrst_file))
simulation.context.setPositions(rst.coordinates[0] * angstrom)
simulation.context.setVelocities(rst.velocities[0] * angstrom / picosecond)
simulation.context.setPeriodicBoxVectors(
    [rst.box[0] * angstrom, 0, 0],
    [0, rst.box[1] * angstrom, 0],
    [0, 0, rst.box[2] * angstrom]
)
simulation.context.setTime(rst.time)
register = StateRegister(simulation)

# Run simulations
for sim in sims:
    simulation.context.reinitialize()
    register.reassign(simulation)

    step = str(sim)
    
    # Set reporters
    outreporter = StateDataReporter(
        os.path.join(path, f'prod_{step}.out'),
        reportInterval=advanced_steps//1000,
        step=True, 
        time=True, 
        potentialEnergy=True, 
        temperature=True,
        kineticEnergy=True, 
        totalEnergy=True, 
        speed=True,
        progress=True, 
        elapsedTime=True, 
        remainingTime=True,
        separator=',', 
        totalSteps=advanced_steps//10,
    )
    dcdreporter = DCDReporter(
        os.path.join(path, f'prod_{step}.dcd'),
        reportInterval=advanced_steps//frames,
    )
    simulation.reporters = [outreporter, dcdreporter]

    # Perform simulation steps
    simulation.step(advanced_steps//10)

    # Save restart file
    rstreporter = RestartReporter(
        os.path.join(path, f'prod_{step}.ncrst'), advanced_steps//10,
        write_multiple=False,
        netcdf=True,
        write_velocities=True
    )
    state = simulation.context.getState(
        getPositions=True, getVelocities=True, enforcePeriodicBox=True
    )
    rstreporter.report(simulation, state)
