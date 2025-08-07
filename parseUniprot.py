import pandas as pd 
import os 
from src.utils.setting import ROOT_DIR
import json
from src.UniProt import searchUniprot
from src.utils.logger import LOGGER
import datetime

####### Set up log file
current_time = datetime.datetime.now().strftime("%Y%m%d-%H%M")
logfile = os.path.join(ROOT_DIR, f'logs/parseUniprot/{current_time}.log')
os.makedirs(os.path.join(ROOT_DIR, f'logs/parseUniprot'), exist_ok=True)
LOGGER.start(logfile)
#######

####### Import Uniprot 
uniprotid = os.path.join(ROOT_DIR, 'data/uniprotid.csv')
df = pd.read_csv(uniprotid)
id_list = df.id.to_list()[:10]
#######

####### Parsing
data = []
for id in id_list:
    try:
        u = searchUniprot(id)
        _dict = {
            'id': u.getAccession(),
            'name': u.getName(),
            'protein': u.getProtein(),
            'gene': u.getGene(),
            'organism': u.getOrganism(),
            'sequence': u.getSequence(),
            'cell_location': u.getCellLocation(),
            'cofactor': u.getCofactor(),
            'binding_site': u.getBindingSite(),
            'active_site': u.getActivateSite(),
            'dna_binding': u.getDNAbinding(),
            'zinc_finger': u.getZincFinger(),
            'pdb': u.getPDBs(),
            'alphafold': u.getAlphaFold(),
        }
        data.append(_dict)
    except Exception as e:
        LOGGER.warn(f'Error while parsing {id}: {e}')
        continue
#######

for label in LOGGER._reports:
    LOGGER.info(f"  {label}: {LOGGER._reports[label]:.2f}s ({LOGGER._report_times[label]} time(s))")
LOGGER.report('Run time elapsed in %.2fs.', "_runtime")
LOGGER.close(logfile)

json_string = json.dumps(data)
filename = 'parseUniprot.json'
filepath = os.path.join(data, filename)
with open(filepath, 'w') as f:
    json.dump(data, f)
    

    