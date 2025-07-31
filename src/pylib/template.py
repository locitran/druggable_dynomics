tleap_template = '''\
set default PBRadii mbondi3
source leaprc.protein.ff14SB
source leaprc.gaff2
loadAmberParams GLYCAM_06j.dat
source leaprc.water.tip3p
loadOff {protein_lib}
solvatebox protein TIP3PBOX 10 iso 
addIonsRand protein Na+ {Na_num_total}
addIonsRand protein Cl- {Cl_num_total}
saveamberparm protein model.prmtop model.inpcrd
savepdb protein model.pdb
quit 
'''.format


mmpbsa_bash_template = '''\
source {AMBERHOME}/amber.sh
{ANTEMMPBSA_EXE} -p {top} -c cmp.prmtop -s :Na+,Cl-,WAT 
mpirun --allow-run-as-root -np 12 {MMPBSAMPI_EXE} -O -i mmpbsa.in -o mmpbsa.out -sp {top} -cp cmp.prmtop -y {traj} -do FINAL_DECOMP_MMPBSA.dat > mmpbsa_mpi.log 2>&1
python {SRCDIR}/pylib/MMPBSA_data_collector.py
'''.format

mmgbsa_template = '''\
&general
   startframe={startframe}, endframe={endframe},
   netcdf=1, verbose=2, keep_files=1,
   entropy=0
/

&gb
   igb=8,
   saltcon={saltcon},
/

&decomp
   csv_format=1, dec_verbose=1, idecomp=2
/
'''.format

mmpbsa_template = '''\
&general
   startframe={startframe}, endframe={endframe},
   netcdf=1, verbose=2, keep_files=1,
   entropy=0
/

&pb
   istrng={saltcon},
   radiopt=0
/

&decomp
   csv_format=1, dec_verbose=1, idecomp=2
/
'''.format

nmode = '''\
Input file for running PB and GB
&general
  interval=1,
  keep_files=0,
  netcdf=1,
/
&nmode
  nmstartframe={startframe},
  nmendframe={endframe},
  nmode_igb=0,
/'''.format

qsub_src = '''#PBS -N Rose{resid}
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1 
#PBS -j oe
#PBS -o job_stdout_filename
cd $PBS_O_WORKDIR

{cartesian_ddg_GCC} @flags_monomer -s {input_pdb} > out
'''.format

sbatch_src = '''#!/bin/bash
#SBATCH --job-name=Rose{resid}
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=example@example.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --output=slurm_openMM%j.log
#SBATCH --partition=COMPUTE1Q
#SBATCH --account=yanglab
#SBATCH --time=0-6

{cartesian_ddg_GCC} @flags_monomer -s {input_pdb} 
'''.format
