import os
import sys
import getopt
import pathlib
import numpy as np
import MDAnalysis as mda
from MDAnalysis.core import groups

CURRENT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.append(str(CURRENT_DIR/'pylib'))

from MDA_PCA import SuperimposePCA ,do_pca
from silhouette import optimal_n_clustering, do_clustering_base_on_silhouette_score

if __name__ == "__main__":

    #input: top and traj file, name
    #output: 2 graph (silhouette score, clustering result base on silhouette score), pdb files of each cluster center

    def usage():
        print('usage: python do_PCA_and_clustering.py --ps "selection" -n output_name --traj traj_1.dcd,traj_2.dcd,traj_3.dcd --top top.prmtop --os "selection" --xray path/to/xray')
        print('\t--ps:\tatom selection in MDAnalysis for doing PCA such as "protein and not (name H*)"')
        print('\t--os:\tatom selection in MDAnalysis for writting output Please notice that you CANNOT appoint -s with select the atom not included in --ps selection.')
        print('\t-s:\toutput cluster center  structures are superimposed. Defulat not.')
        print('\t-n:\toutput file name, defalut: protein')
        print('\t--traj:\ttrajectory file each trajectory file split with "," and without space')
        print('\t--top:\ttopology file in MDAnalysis readable format (prmtop, parm7...)')
        print('\t--xray:\tpath/to/xray_pdb: add xray structure in clustering')

    opt_list, args = getopt.getopt(sys.argv[1:], 'n:hs' ,["traj=","top=","xray=","ps=","os="])
    opt_list = dict(opt_list)
    
    if '-h' in opt_list:
        usage()
        sys.exit(2)

    if '-s' in opt_list:
        superimpose_flag = True
    else:
        superimpose_flag = False
  
    opt_list = dict(opt_list)
    
    trajectory_files = []
    for traj in opt_list['--traj'].split(','):
        trajectory_files.append(os.path.abspath(traj))

    topology_file = os.path.abspath(opt_list['--top'])
    PCA_selection = opt_list['--ps']

    if '-n' in opt_list:
        filename = os.path.join(os.path.dirname(topology_file), opt_list['-n'])
    else:
        filename = os.path.join(os.path.dirname(topology_file), 'protein')

    if '--xray' in opt_list:
        xray = os.path.abspath(opt_list.get('--xray'))
        universe = mda.Universe(topology_file, xray, *trajectory_files)
        xray_flag = True
    else:
        universe = mda.Universe(topology_file, *trajectory_files)
        xray_flag = False

    if '--os' in opt_list:
        output_selection = opt_list["--os"]
    else:
        output_selection = PCA_selection

    pc1, pc2, pc1_var, pc2_var, _ , _, new_positions = do_pca(universe,PCA_selection)
    print (pc1, pc2)
    pc1, pc2 = np.real(pc1), np.real(pc2)
    do_clustering_base_on_silhouette_score(universe, [pc1, pc2], [pc1_var, pc2_var], filename, xray_flag, output_selection, new_positions, superimpose_flag)
    #print(universe, [pc1, pc2], [pc1_var, pc2_var], filename, xray_flag, output_selection, new_positions, superimpose_flag)