import sys
sys.path.append('${{AMBERPYTHIONLIB}}')
sys.path.append('${{PARMEDPYTHONLIB}}')

import numpy as np
import json
from MMPBSA_mods import API as MMPBSA_API

# MMPBSAAnalyzer
def gb_diff(data):
    '''
    Differences (Complex - Receptor - Ligand),
        Average DELTA TOTAL
    '''
    c = data['gb']['complex']['TOTAL']
    r = data['gb']['receptor']['TOTAL']
    l = data['gb']['ligand']['TOTAL']

    #return np.mean(c - r - l)
    return c - r - l


def pb_diff(data):
    '''
    Differences (Complex - Receptor - Ligand),
        Average DELTA TOTAL
    '''
    c = data['pb']['complex']['TOTAL']
    r = data['pb']['receptor']['TOTAL']
    l = data['pb']['ligand']['TOTAL']

    #return np.mean(c - r - l)
    return c - r - l


mmpbsa_data = MMPBSA_API.load_mmpbsa_info('_MMPBSA_info')
data = dict()
data['gb'] = gb_diff(mmpbsa_data).tolist()
#data['pb'] = pb_diff(mmpbsa_data).tolist()

json_file = 'mmpbsa_info.json'
with open(json_file, 'w') as f:
    f.write(json.dumps(data))
    
# /usr/bin/python2.7 mmpbsa_info.py

