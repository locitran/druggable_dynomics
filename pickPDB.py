import os 
import json
import datetime

from src.utils.setting import ROOT_DIR
from src.UniProt import searchUniprot, pickPDBfromUniprot
from src.utils.logger import LOGGER

####### Set up log file
current_time = datetime.datetime.now().strftime("%Y%m%d-%H%M")
logfile = os.path.join(ROOT_DIR, f'logs/pickPDB/{current_time}.log')
os.makedirs(os.path.join(ROOT_DIR, f'logs/pickPDB'), exist_ok=True)
LOGGER.start(logfile)
#######

####### Import Uniprot 
# uniprotid = os.path.join(ROOT_DIR, 'data/uniprotid.csv')
# df = pd.read_csv(uniprotid)
# id_list = df.id.to_list()
idlist = ['P05106','P05093','O14764','P18505','P08588','P23219','P19793','P35367','P31645','A8TX70','P02462']
#######

####### Parsing
data = []
for id in idlist:
    try:
        u = searchUniprot(id)
        data = {
            'id': u.getAccession(),
            'name': u.getName(),
            'protein': u.getProtein(),
            'gene': u.getGene(),
            'organism': u.getOrganism(),
            'sequence': u.getSequence(),
            'cell_location': u.getCellLocation(),
            'cofactor': u.getCofactor(),
            'binding_site': u.getBindingSite(),
            'active_site': u.getActiveSite(),
            'dna_binding': u.getDNAbinding(),
            'zinc_finger': u.getZincFinger(),
            'site': u.getSite(),
            'pdb': u.getPDBs(),
            'alphafold': u.getAlphaFold(),
            }
        PDBrank = pickPDBfromUniprot(data, to_file=f'data/{id}_PDBrank.csv')
    except Exception as e:
        LOGGER.warn(f'Error while parsing {id}: {e}')
        continue
#######

for label in LOGGER._reports:
    LOGGER.info(f"  {label}: {LOGGER._reports[label]:.2f}s ({LOGGER._report_times[label]} time(s))")
LOGGER.report('Run time elapsed in %.2fs.', "_runtime")
LOGGER.close(logfile)

# json_string = json.dumps(data)
# filename = 'parseUniprot_site_881.json'
# filepath = os.path.join(ROOT_DIR, 'data', filename)
# with open(filepath, 'w') as f:
#     json.dump(data, f)
    

    