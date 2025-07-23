import os
import MDAnalysis as mda
from openmm import *
from openmm.app import *
from openmm.unit import *
from parmed.openmm import NetCDFReporter, RestartReporter
from OpenMMUtilities import *

def do_openMM_simulation(model_folder, prod_time, frames, GPU_id):

    prmfile = os.path.join(model_folder, 'model.prmtop')
    crdfile = os.path.join(model_folder, 'model.inpcrd')
    prmtop = AmberPrmtopFile(prmfile)
    inpcrd = AmberInpcrdFile(crdfile)
    prod_time = prod_time * nanoseconds

    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=10*angstroms,
        constraints=HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )

    dt = 1.0 * femtoseconds
    fric_coeff = 1.0 / picosecond #frictionCoeff
    temperature = 0.0
    integrator = LangevinIntegrator(temperature, fric_coeff, dt)
    #integrator.setConstraintTolerance(0.00001)

    # platform settings
    platform = Platform.getPlatformByName('CUDA')
    platformProp = dict([
        ('CudaDeviceIndex', GPU_id),
        ('CudaPrecision', 'mixed'),
    ])

    # initialize simulation
    simulation = Simulation(prmtop.topology, system, integrator, platform, platformProp)
    simulation.context.setPositions(inpcrd.positions)
    #simulation.context.setVelocitiesToTemperature(0.0*kelvin)
    simulation.context.setPeriodicBoxVectors(*inpcrd.getBoxVectors()) 

    # initialze register
    register = StateRegister(simulation)

    u = mda.Universe(prmfile,crdfile)
    protein_Universe = u.select_atoms('protein and not name H*')
    backbone_Universe = u.select_atoms('backbone')
    sidechain_Universe = u.select_atoms('protein and not name C CA O N H*')
    
    restraintmasks = dict([

        ('Protein', protein_Universe.atoms.indices), #@* & !(@H=|:WAT@O|@Na+|@Cl-)
        
        ('Backbone', backbone_Universe.atoms.indices), 
        
        ('Sidechain', sidechain_Universe.atoms.indices),

    ])

    # Minimization_1
    '''
    Restraint force constant: 9.56 kcal/mol/A**2
    Restraint residues heavy atoms
    Minimize waters, ions, hydrogens
    500 steps
    '''
    restraint_wt = 9.56 * kilocalorie/mole/angstrom ** 2
    restraintmask = restraintmasks['Protein']
    restraint_ref = inpcrd.positions

    restraintForce = HarmonicRestraint(restraint_wt, restraintmask, restraint_ref)
    restraintForceIndex = system.addForce(restraintForce)

    simulation.context.reinitialize()
    register.reassign(simulation)

    simulation.minimizeEnergy(maxIterations=500)

    register.reserve(simulation)
    system.removeForce(restraintForceIndex)

    # Minimization_2
    '''
    Restraint force constant: 0 kcal/mol/A**2
    Minimize residues side chains
    500 steps
    '''

    simulation.context.reinitialize()
    register.reassign(simulation)

    simulation.minimizeEnergy(maxIterations=500)

    register.reserve(simulation)

    # Heating
    '''
    Restraint force constant: 1 kcal/mol/A**2
    Restraint residues CA atoms
    NVT ensemble
    Heat to 320.0 K
    Start from 100.0 K
    1000 picoseconds (500000 steps * 0.002 fs/step)
    '''
    Heating_320K_time = 0.5 * nanoseconds
    Heating_320K_steps = Heating_320K_time / dt
    Heating_310K_time = 2.125 * nanoseconds
    gradient_time = 1.0 * picoseconds
    gradient_steps = gradient_time / dt
    Heating_310K_steps = Heating_310K_time / dt    #for protein comformation exploring 
    tempi = 10.15 * kelvin
    restraint_wt_bb = 4.78 * kilocalorie/mole/angstrom ** 2
    restraintmask_bb = restraintmasks['Backbone']
    restraint_wt_sc = 2.39 * kilocalorie/mole/angstrom ** 2
    restraintmask_sc = restraintmasks['Sidechain']   
    restraint_ref = inpcrd.positions

    restraintForce_bb = HarmonicRestraint(restraint_wt_bb, restraintmask_bb, restraint_ref)
    restraintForceIndex_bb = system.addForce(restraintForce_bb)
    restraintForce_sc = HarmonicRestraint(restraint_wt_sc, restraintmask_sc, restraint_ref)
    restraintForceIndex_sc = system.addForce(restraintForce_sc)

    simulation.context.setVelocitiesToTemperature(tempi)
    simulation.context.getIntegrator().setTemperature(tempi)

    simulation.context.reinitialize()
    register.reassign(simulation)

    # set reporter
    outreporter = StateDataReporter(
         os.path.join(model_folder, 'heat.out'),
        reportInterval=5000,
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
        totalSteps=Heating_320K_steps + Heating_310K_steps,
    )
    ncreporter = DCDReporter(
         os.path.join(model_folder, 'heat.dcd'),
        reportInterval=5000,
    )
    simulation.reporters = [outreporter, ncreporter]
    simulation.currentStep = 0 # let simulation step back to zero

    for grad in range(11):
        simulation.step(int(gradient_steps))
        if grad == 0 or grad == 10:   
            tempi += 20 * kelvin
        else:
            tempi += 30 * kelvin
        simulation.context.getIntegrator().setTemperature(tempi)


    Heating_320K_steps -= int(gradient_steps) * 11
    simulation.integrator.setTemperature(320.15 * kelvin)
    simulation.step(Heating_320K_steps)

    simulation.integrator.setTemperature(310.15 * kelvin)
    simulation.step(Heating_310K_steps)

    register.reserve(simulation)
    system.removeForce(restraintForceIndex_sc)
    system.removeForce(restraintForceIndex_bb)

    # Equilibration_NPT_1
    '''
    Restraint force constant: 2.39 kcal/mol/A**2
    Restraint residues backbone atoms
    Restraint force constant: 1.20 kcal/mol/A**2
    Restraint residues sidechain atoms
    NPT ensemble
    Maintain 310.15 K
    Maintain 1 atmospheres
    125 picoseconds 
    '''
    equilibration_NPT_time = 0.125 * nanoseconds
    equilibration_NPT_step = equilibration_NPT_time / dt

    temp = 310.15 * kelvin
    pressure = 1 * atmospheres
    restraint_wt_bb = 2.39 * kilocalorie/mole/angstrom ** 2
    restraintmask_bb = restraintmasks['Backbone']
    restraint_wt_sc = 1.20 * kilocalorie/mole/angstrom ** 2
    restraintmask_sc = restraintmasks['Sidechain']   
    restraint_ref = inpcrd.positions

    simulation.context.getIntegrator().setTemperature(temp)
    system.addForce(MonteCarloBarostat(pressure, temp, 25)) # frequency=25

    restraintForce_bb = HarmonicRestraint(restraint_wt_bb, restraintmask_bb, restraint_ref)
    restraintForceIndex_bb = system.addForce(restraintForce_bb)
    restraintForce_sc = HarmonicRestraint(restraint_wt_sc, restraintmask_sc, restraint_ref)
    restraintForceIndex_sc = system.addForce(restraintForce_sc)

    simulation.context.reinitialize()
    register.reassign(simulation)

    # set reporter
    outreporter = StateDataReporter(
         os.path.join(model_folder, 'equi.out'),
        reportInterval=7500,
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
        totalSteps=750000,
    )
    ncreporter = DCDReporter(
         os.path.join(model_folder, 'equi.dcd'),
        reportInterval=7500,
    )
    simulation.reporters = [outreporter, ncreporter]
    simulation.currentStep = 0 # let simulation step back to zero

    simulation.step(equilibration_NPT_step)

    register.reserve(simulation)
    system.removeForce(restraintForceIndex_sc)
    system.removeForce(restraintForceIndex_bb)

    # Equilibration_NPT_2
    '''
    Restraint force constant: 1.20 kcal/mol/A**2
    Restraint residues backbone atoms
    Restraint force constant: 0.48 kcal/mol/A**2
    Restraint residues sidechain atoms
    NPT ensemble
    Maintain 310.15 K
    Maintain 1 atmospheres
    250 picoseconds 
    '''

    dt = 2.0 * femtoseconds # change time step to 2
    # fric_coeff = 1.0 / picosecond #frictionCoeff
    temperature = 310.0 * kelvin
    # integrator = LangevinIntegrator(temperature, fric_coeff, dt)

    # simulation = Simulation(prmtop.topology, system, integrator, platform, platformProp)
    # simulation.context.reinitialize()
    # register.reassign(simulation)

    equilibration_NPT_time = 0.250 * nanoseconds
    equilibration_NPT_step = equilibration_NPT_time / dt

    temp = 310.15 * kelvin
    pressure = 1 * atmospheres
    restraint_wt_bb = 1.20 * kilocalorie/mole/angstrom ** 2
    restraintmask_bb = restraintmasks['Backbone']
    restraint_wt_sc = 0.48 * kilocalorie/mole/angstrom ** 2
    restraintmask_sc = restraintmasks['Sidechain']   
    restraint_ref = inpcrd.positions

    restraintForce_bb = HarmonicRestraint(restraint_wt_bb, restraintmask_bb, restraint_ref)
    restraintForceIndex_bb = system.addForce(restraintForce_bb)
    restraintForce_sc = HarmonicRestraint(restraint_wt_sc, restraintmask_sc, restraint_ref)
    restraintForceIndex_sc = system.addForce(restraintForce_sc)

    simulation.context.getIntegrator().setStepSize(dt)
    # simulation.context.getIntegrator().setTemperature(temp)
    # system.addForce(MonteCarloBarostat(pressure, temp, 25)) # frequency=25

    simulation.context.reinitialize()
    register.reassign(simulation)

    simulation.step(equilibration_NPT_step)

    register.reserve(simulation)
    system.removeForce(restraintForceIndex_sc)
    system.removeForce(restraintForceIndex_bb)

    # Equilibration_NPT_3
    '''
    Restraint force constant: 0.48 kcal/mol/A**2
    Restraint residues backbone atoms
    Restraint force constant: 0.12 kcal/mol/A**2
    Restraint residues sidechain atoms
    NPT ensemble
    Maintain 310.15 K
    Maintain 1 atmospheres
    500 picoseconds 
    '''
    equilibration_NPT_time = 0.500 * nanoseconds
    equilibration_NPT_step = equilibration_NPT_time / dt

    temp = 310.15 * kelvin
    pressure = 1 * atmospheres
    restraint_wt_bb = 0.48 * kilocalorie/mole/angstrom ** 2
    restraintmask_bb = restraintmasks['Backbone']
    restraint_wt_sc = 0.12 * kilocalorie/mole/angstrom ** 2
    restraintmask_sc = restraintmasks['Sidechain']   
    restraint_ref = inpcrd.positions

    restraintForce_bb = HarmonicRestraint(restraint_wt_bb, restraintmask_bb, restraint_ref)
    restraintForceIndex_bb = system.addForce(restraintForce_bb)
    restraintForce_sc = HarmonicRestraint(restraint_wt_sc, restraintmask_sc, restraint_ref)
    restraintForceIndex_sc = system.addForce(restraintForce_sc)

    # simulation.context.getIntegrator().setTemperature(temp)
    # system.addForce(MonteCarloBarostat(pressure, temp, 25)) # frequency=25

    simulation.context.reinitialize()
    register.reassign(simulation)

    simulation.step(equilibration_NPT_step)

    register.reserve(simulation)
    system.removeForce(restraintForceIndex_sc)
    system.removeForce(restraintForceIndex_bb)

    # Equilibration_NPT_4
    '''
    Restraint force constant: 0.12 kcal/mol/A**2
    Restraint residues backbone atoms
    NPT ensemble
    Maintain 310.15 K
    Maintain 1 atmospheres
    500 picoseconds 
    '''
    equilibration_NPT_time = 0.500 * nanoseconds
    equilibration_NPT_step = equilibration_NPT_time / dt

    temp = 310.15 * kelvin
    pressure = 1 * atmospheres
    restraint_wt_bb = 0.12 * kilocalorie/mole/angstrom ** 2
    restraintmask_bb = restraintmasks['Backbone']
    restraint_ref = inpcrd.positions

    restraintForce_bb = HarmonicRestraint(restraint_wt_bb, restraintmask_bb, restraint_ref)
    restraintForceIndex_bb = system.addForce(restraintForce_bb)

    # simulation.context.getIntegrator().setTemperature(temp)
    # system.addForce(MonteCarloBarostat(pressure, temp, 25)) # frequency=25

    simulation.context.reinitialize()
    register.reassign(simulation)

    simulation.step(equilibration_NPT_step)

    register.reserve(simulation)
    system.removeForce(restraintForceIndex_bb)

    # Production
    '''
    NPT ensemble
    Maintain 310.0 K
    Maintain 1 atmospheres
    prod_time = 10 nanoseconds (5,000,000 steps * 2 fs(dt)/step = 10,000,000 fs = 10,000 ps)
    snapshot per 100 ps (50,000 steps)
    '''
    # temp = 310.0 * kelvin
    # pressure = 1 * atmospheres

    # simulation.context.getIntegrator().setTemperature(temp)
    # system.addForce(MonteCarloBarostat(pressure, temp, 25))

    simulation.currentStep = 0 # let simulation step back to zero
    simulation.context.setTime(0.0) # let simulation time back to zero

    advanced_steps = int(prod_time / dt)
    # 5,000,000 steps = 1000,000,000 fs / 2 fs
    
    for part in range(1, 11):

        simulation.context.reinitialize()
        register.reassign(simulation)

        # set reporter
        outreporter = StateDataReporter(
            os.path.join(model_folder, f'prod_{part}.out'),
            reportInterval=advanced_steps//1000, # 5,000,000 steps / 1000 = 5000 steps *2fs = 10,000 fs = 10 ps
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
        ncreporter = DCDReporter(
            os.path.join(model_folder, f'prod_{part}.dcd'),
            reportInterval=advanced_steps//frames, # 5,000,000 steps / 2000 = 5000 steps *2fs = 10,000 fs = 10 ps
        )
        simulation.reporters = [outreporter, ncreporter]

        simulation.step(advanced_steps//10)
            
        register.reserve(simulation)
        #system.removeForce(restraintForceIndex)

        # save restart file
        rstreporter = RestartReporter(
            os.path.join(model_folder, f'prod_{part}.ncrst'), advanced_steps//10,
            write_multiple=False,
            netcdf=True,
            write_velocities=True
        )
        state = simulation.context.getState(
            getPositions=True, getVelocities=True, enforcePeriodicBox=True
        )
        rstreporter.report(simulation, state)
