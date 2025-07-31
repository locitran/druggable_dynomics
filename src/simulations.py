import re
import os
import sys
import getopt
import dotenv 
import pathlib
import subprocess as sp

CURRENT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.append(str(CURRENT_DIR/'pylib'))

from template import *
from simulation_OpenMM import do_openMM_simulation
from cal_required_ions_number_for_concentration import count_ion
from get_leaplog_water_ion_numbers import get_leaplog_water_ion_numbers

def do_simulation(work_folder, protein_folder, ion_concentration, prod_time, frames, GPU_id):
    
    tleapExe = os.environ['TLEAP_EXE']

    md_folder = os.path.join(work_folder, 'MD')
    log_folder = os.path.join(md_folder, 'Log')

    os.makedirs(md_folder,exist_ok=True)
    os.makedirs(log_folder,exist_ok=True)

    protein_lib = os.path.join(protein_folder,'protein.lib')

    script = tleap_template(protein_lib=protein_lib,
                            Na_num_total=0,
                            Cl_num_total=0)

    with open(os.path.join(md_folder,'tleap.in'),'w') as f:
        f.write(script)

    command = f'{tleapExe} -f {md_folder}/tleap.in'
    child = sp.Popen(command.split(), stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, cwd=md_folder)
    stdout, stderr = child.communicate()
    child.terminate()

    if stderr:
        with open(os.path.join(log_folder,'tleap.err'),'w') as f:
            f.write(stderr.decode())
    
    water_num, Na_num, Cl_num = get_leaplog_water_ion_numbers(md_folder + '/leap.log')
    N_total_Na, N_total_Cl = count_ion(water_num,ion_concentration , Na_num, Cl_num)

    script = tleap_template(protein_lib=protein_lib,
                            Na_num_total=N_total_Na,
                            Cl_num_total=N_total_Cl)

    with open(os.path.join(md_folder,'tleap.in'),'w') as f:
        f.write(script)

    command = f'{tleapExe} -f {md_folder}/tleap.in'
    child = sp.Popen(command.split(), stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, cwd=md_folder)
    stdout, stderr = child.communicate()
    child.terminate()

    if stderr:
        with open(os.path.join(log_folder,'tleap.err'),'a') as f:
            f.write(stderr.decode())
        
    do_openMM_simulation(md_folder, prod_time, frames, GPU_id)

if __name__ == '__main__':
    
    import time
    start_time = time.time()
    
    def usage():
        print('usage: python simulations.py -o output_folder -g gpu_id -i protein_folder -c ion_concentration -t production_time -f frames')
        print('\t-o:\tworking folder where all the output files will be located.')
        print('\t-i:\tthe folder which is produced by prepare_protein.py')
        print('\t-c:\tthe NaCl ion concentration (mol/L)')
        print('\t-t:\tthe production time in ns, 1000 ns')
        print('\t-f:\tthe total frames you will get after MD simulation, 2000')
        print('\t-g:\tgpu id')
        print('\t-h:\tprint help message')

    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'o:i:c:g:t:f:h')
        opt_list = dict(opt_list)
        if '-h' in opt_list:
            usage()
            sys.exit(2)

        protein_folder = os.path.abspath(opt_list['-i'])
        ion_concentration = opt_list['-c']
        GPU_id = opt_list['-g']
        prod_time = int(opt_list['-t'])
        frames = int(opt_list['-f'])

    except Exception as e:
        print(e)
        usage()
        sys.exit(2)

    if '-o' in opt_list:
        output_folder = os.path.abspath(opt_list['-o'])       
    else:
        output_folder = os.path.dirname(protein_folder)

    CURRENT_DIR = pathlib.Path(__file__).resolve().parent
    ENV_FILE_PATH = CURRENT_DIR / '.env'
    dotenv.load_dotenv(str(ENV_FILE_PATH))

    do_simulation(output_folder, protein_folder, ion_concentration, prod_time, frames, GPU_id)
    print("--- %.2f seconds ---" % (time.time() - start_time))
