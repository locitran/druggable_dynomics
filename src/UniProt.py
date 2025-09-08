"""This module defines a class and relative functions for parsing Uniprot
information and corresponding PDB.

Some important definitions should be acknowledged. 

Reference: https://www.rcsb.org/docs/general-help/identifiers-in-pdb
1. PDB entity: (Entity level Identifier)
An entity is a chemically distinct part of a structure. One entry may have multiple copies (instances) of a given entity

2. PDB instance: (Instance level Identifiers)
An instance is a distinct copy of an entity or molecule.

For example, a heterotetramer hemoglobin (4hhb) contains two copies of the hemoglobin alpha chain (or chains A and C) 
and two copies of the beta chain (or chains B and D). People then name it as follow:
- the hemoglobin alpha: entity 1 --> chains A and C: instances A and C
- the hemoglobin beta: entity 2 --> chains B and D: instances B and D
"""

import time
import re

from prody import queryUniprot
from prody.utilities import openURL
from rcsbapi.data import DataQuery as Query
from rcsbapi.search import AttributeQuery, NestedAttributeQuery

from src.utils.logger import LOGGER

comma_splitter = re.compile(r'\s*,\s*').split
ns = {'up': 'http://uniprot.org/uniprot'}

COFACTOR = {
    'CHEBI:17996': ['CL'], # chloride
    'CHEBI:18420': ['MG'], # Mg(2+)
    'CHEBI:25516': ['NI'], # Ni cation
    'CHEBI:29033': ['FE2'],# Fe(2+)
    'CHEBI:24875': ['FE2', 'FE'], # Fe cation
    'CHEBI:29035': ['MN'], # Mn(2+)
    'CHEBI:29036': ['CU'], # Cu(2+)
    'CHEBI:29101': ['NA'],  # Na(+)
    'CHEBI:29103': ['K'],  # K(+)
    'CHEBI:29105': ['ZN'],  # Zn(2+)
    'CHEBI:29108': ['CA'],  # Ca(2+)
    'CHEBI:60240': ['FE2', 'MN', 'CA', 'CU', 'MG', 'NI', 'ZN'], # a divalent metal cation
    'CHEBI:30413': ['HEM'], # heme
    'CHEBI:60344': ['COH', 'HEB'], # heme b
    'CHEBI:190135':['FES'], # [2Fe-2S] cluster
    'CHEBI:49883': ['SF4'],     # [4Fe-4S] cluster
    'CHEBI:71302': ['MOM', 'MTE'], # Mo-molybdopterin
    'CHEBI:57692': ['FAD'], # FAD
    'CHEBI:58210': ['FMN'], # FMN
    'CHEBI:58937': ['TPP'], # thiamine diphosphate
    'CHEBI:597326': ['PLP'], # pyridoxal 5'-phosphate
    'CHEBI:unk_0': ['NAD', 'NDP', 'NAP'], # NICOTINAMIDE-ADENINE-DINUCLEOTIDE
    'CHEBI:unk_1': ['COA'], # COENZYME A
    'CHEBI:unk_2': ['B1Z', 'B12'], # Adenosylcobalamin or COBALAMIN # 5'-Deoxyadenosylcobalamin (vitamin B12)
    'CHEBI:unk_3': ['BTN'], # Biotin (biocytin)
    'CHEBI:unk_4': ['THG'], # Tetrahydrofolate (THF)
}

COFACTOR_COENZYME = ['FAD', 'FMN', 'TPP', 'PLP', 'NAD', 'NDP', 'NAP', 'COA', 'B1Z', 'B12', 'BTN', 'THG']

class UniprotRecord(object):
    """Wrapper for UniProt data with methods for accessing fields and parsing related PDB entries."""
    
    def __init__(self, data):
        self._rawdata = data
        self._pdbdata = []
        self._parse()

    def __repr__(self):
        return '<UniprotRecord: %s>' % self.getTitle()

    def __str__(self):
        return self.getTitle()

    def setData(self, value):
        self._rawdata = value
        self._parse()

    def getData(self):
        return self._rawdata

    def getPDBs(self):
        return self._pdbdata

    def getAccession(self, index=0):
        """Get primary accession."""
        return self.getEntry('accession', index)

    def getName(self, index=0):
        """Get entry name."""
        return self.getEntry('name', index)

    def getProtein(self, index=0):
        """Get recommended and alternative protein names."""
        protein = self.getEntry('protein', index)
        try:
            recommend_elem = protein.find('up:recommendedName', ns)
            alternative_elem = protein.find('up:submittedName', ns)
            
            recommend_name = recommend_elem.find('up:fullName', ns) if recommend_elem is not None else None
            alter_fullname = alternative_elem.find('up:fullName', ns) if alternative_elem is not None else None
            alter_shortname = alternative_elem.find('up:shortName', ns) if alternative_elem is not None else None
            
            recommend_name = recommend_name.text if recommend_name is not None else None
            alter_fullname = alter_fullname.text if alter_fullname is not None else None
            alter_shortname = alter_shortname.text if alter_shortname is not None else None
        except:
            # Fallback if structure differs
            submitted_name = protein.find('up:submittedName/up:fullName', ns)
            recommend_name = submitted_name.text if submitted_name is not None else None
            alter_fullname = None
            alter_shortname = None
        
        return {
            'recommend_name': recommend_name,
            'alter_fullname': alter_fullname,
            'alter_shortname': alter_shortname            
        }

    def getGene(self, index=0):
        """Get primary gene name."""
        try:
            gene = self.getEntry('gene', index)
            name_elem = gene.find('up:name[@type="primary"]', ns)
            return name_elem.text if name_elem is not None else None
        except Exception as e:
            LOGGER.warn(f'Error parsing gene name: {e}')
            return None

    def getOrganism(self, index=0):
        """Get organism scientific and common names, taxonomy."""
        organism = self.getEntry('organism', index)
        sci_name = organism.find('up:name[@type="scientific"]', ns)
        com_name = organism.find('up:name[@type="common"]', ns)
        db_ref = organism.find('up:dbReference[@type="NCBI Taxonomy"]', ns)
        lineage_tags = organism.findall('up:lineage/up:taxon', ns)
        return {
            'scientific_name': sci_name.text.strip() if sci_name is not None else None,
            'common_name': com_name.text.strip() if com_name is not None else None,
            'taxonomy_id': db_ref.attrib['id'] if db_ref is not None else None,
            'lineage': [taxon.text.strip() for taxon in lineage_tags if taxon.text]
        }

    def getCellLocation(self):
        return self._cell_location

    def getDNAbinding(self):
        return self._dna_binding

    def getZincFinger(self):
        return self._zinc_finger

    def getActiveSite(self):
        return self._active_site

    def getBindingSite(self):
        return self._binding_site

    def getSite(self):
        return self._site

    def getAlphaFold(self):
        """Return AlphaFold DB ID if present."""
        AlphaFoldDB = None
        for key, value in self._rawdata.items():
            if not key.startswith('dbReference'):
                continue
            if type(value) != list or len(value) != 2:
                continue
            # List like: [('type','AlphaFoldDB'), ('id','QXXXXX')]
            if value[0][1] == 'AlphaFoldDB':
                AlphaFoldDB = value[1][1]
                break
        return AlphaFoldDB

    def getCofactor(self):
        return self._cofactors

    # New getters for additional parsed features:
    def getPTMProcessing(self):
        return self._ptm_processing

    def getInteractions(self):
        return self._interactions

    def getFamilyDomains(self):
        return self._family_domains

    def getComplexPortals(self):
        """Return list of complex names from ComplexPortal references."""
        return self._complex_portals
            
    def getTitle(self):
        uid = self.getAccession()
        name = self.getName()
        return f'{uid} ({name})'

    def getEntry(self, item, index=0):
        key = f'{item}{index:4d}'
        if key in self._rawdata:
            return self._rawdata[key]
        else:
            raise KeyError(f'{item} does not exist in the record')

    def _parseDNAbinding(self):
        data = self._rawdata
        dna_binding = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            if value.get('type') != "DNA-binding region":
                continue
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            dna_binding.append({'description': descp, 'begin': begin, 'end': end})
        self._dna_binding = dna_binding

    def _parseZincfinger(self):
        data = self._rawdata
        zinc_finger = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            if value.get('type') != "zinc finger region":
                continue
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            zinc_finger.append({'description': descp, 'begin': begin, 'end': end})
        self._zinc_finger = zinc_finger

    def _parseActiveSite(self):
        data = self._rawdata
        active_site = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            if value.get('type') != "active site":
                continue
            descp = value.get('description')
            pos_elem = value.find('up:location/up:position', ns)
            pos = int(pos_elem.attrib.get('position')) if pos_elem is not None else None
            active_site.append({'description': descp, 'position': pos})
        self._active_site = active_site

    def _parseBindingSite(self):
        data = self._rawdata
        binding_site = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            if value.get('type') != "binding site":
                continue
            descp = value.get('description', None)
            loc_elem = value.find('up:location', ns)
            pos_elem = loc_elem.find('up:position', ns)
            begin_elem = loc_elem.find('up:begin', ns)
            end_elem = loc_elem.find('up:end', ns)
            if pos_elem is None:
                pos = f"{begin_elem.attrib.get('position')}-{end_elem.attrib.get('position')}"
            else:
                pos = pos_elem.attrib.get('position')
            ligand_elem = value.find('up:ligand', ns)
            ligand_name = ligand_elem.find('up:name', ns)
            ligand_name = ligand_name.text if ligand_name is not None else None
            ligand_chebi = ligand_elem.find('up:dbReference[@type=\"ChEBI\"]', ns)
            ligand_chebi = ligand_chebi.attrib['id'] if ligand_chebi is not None else None
            binding_site.append({
                'position': pos,
                'description': descp,
                'name': ligand_name,
                'chebi': ligand_chebi
            })
        self._binding_site = binding_site

    def _parseSite(self):
        data = self._rawdata
        site = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            if value.get('type') != "site":
                continue
            descp = value.get('description')
            loc_elem = value.find('up:location', ns)
            pos_elem = loc_elem.find('up:position', ns)
            begin_elem = loc_elem.find('up:begin', ns)
            end_elem = loc_elem.find('up:end', ns)
            if pos_elem is None:
                pos = f"{begin_elem.attrib.get('position')}-{end_elem.attrib.get('position')}"
            else:
                pos = pos_elem.attrib.get('position')
            site.append({'position': pos, 'description': descp})
        self._site = site

    def _parseCofactor(self):
        data = self._rawdata
        cofactors = []
        for key, value in data.items():
            if not key.startswith('comment'):
                continue
            if type(value) == list:
                continue
            if value.get('type') != "cofactor":
                continue
            cf_elem = value.find('up:cofactor', ns)
            cf_name = cf_elem.find('up:name', ns)
            cf_chebi = cf_elem.find('up:dbReference[@type=\"ChEBI\"]', ns)
            cf_chebi = cf_chebi.attrib['id'] if cf_chebi is not None else None
            cofactors.append({
                'name': cf_name.text,
                'chebi': cf_chebi
            })
        self._cofactors = cofactors

    def _parseCellLocation(self):
        data = self._rawdata
        cell_location = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            loc_type = value.get('type')
            if loc_type not in ['topological domain', 'transmembrane region', 'intramembrane region']:
                continue
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            cell_location.append({
                'type': loc_type,
                'description': descp,
                'begin': begin,
                'end': end
            })
        self._cell_location = cell_location

    def _parsePTMProcessing(self):
        data = self._rawdata
        ptms = []
        types_of_interest = ["signal peptide", "transit peptide", "chain", "propeptide", "peptide",
                             "disulfide bond", "glycosylation site", "modified residue", "lipidation site"]
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            type_ = value.get('type')
            if type_ not in types_of_interest:
                continue
            descp = value.get('description')
            loc_elem = value.find('up:location', ns)
            if loc_elem is None:
                continue
            pos_elem = loc_elem.find('up:position', ns)
            begin_elem = loc_elem.find('up:begin', ns)
            end_elem = loc_elem.find('up:end', ns)
            if pos_elem is not None:
                pos = pos_elem.attrib.get('position')
                ptms.append({'type': type_, 'position': pos, 'description': descp})
            else:
                begin = begin_elem.attrib.get('position') if begin_elem is not None else None
                end = end_elem.attrib.get('position') if end_elem is not None else None
                ptms.append({'type': type_, 'begin': begin, 'end': end, 'description': descp})
        self._ptm_processing = ptms

    def _parseComplexPortals(self):
        """Return the list of ComplexPortal complex names."""
        data = self._rawdata
        self._complex_portals = []

        for key, value in data.items():
            if not key.startswith('dbReference'):
                continue

            if type(value) is list or type(value) is dict:
                continue

            if value.get('type') != 'ComplexPortal':
                continue
            
            cp_id = value.get('id')
            property_elem = value.find('up:property', ns)
            if property_elem is None:
                continue
            name = property_elem.get('value')
            self._complex_portals.append(
                {
                    'id': cp_id,
                    'name': name
                }
            )

    def _parseFamilyDomains(self):
        data = self._rawdata
        family_domains = []
        types_of_interest = [
            "domain", "region of interest", "coiled-coil region", 
            "short sequence motif", "repeat"
        ]
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            type_ = value.get('type')
            if type_ not in types_of_interest:
                continue
            descp = value.get('description')
            loc_elem = value.find('up:location', ns)
            if loc_elem is None:
                continue
            begin_elem = loc_elem.find('up:begin', ns)
            end_elem = loc_elem.find('up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            family_domains.append({'type': type_, 'description': descp, 'begin': begin, 'end': end})
        self._family_domains = family_domains

    def _parseInteractions(self):
        data = self._rawdata
        interactions = []
        self_acc = self.getAccession()
        for key, value in data.items():
            if not key.startswith('comment'):
                continue
            if value.get('type') != "interaction":
                continue
            interactants = value.findall('up:interactant', ns)
            for ia in interactants:
                ia_id_elem = ia.find('up:id', ns)
                if ia_id_elem is None:
                    continue
                ia_id = ia_id_elem.text
                label_elem = ia.find('up:label', ns)
                ia_label = label_elem.text if label_elem is not None else None
                if ia_id != self_acc:
                    interactions.append({'id': ia_id, 'label': ia_label})
        self._interactions = interactions

    def _parseSeqAnnotfromPDB(self, pdb_instances):
        """pdb_instances ['1AID.A', '1AID.B']"""
        uniprotid = self.getAccession()
        query = Query(
            input_type="polymer_entity_instances",
            input_ids=pdb_instances,
            return_data_list=[
                "polymer_entity_instances.rcsb_id",
                "rcsb_polymer_instance_feature.type",
                "rcsb_polymer_instance_feature.feature_positions.beg_seq_id",
                "rcsb_polymer_instance_feature.feature_positions.end_seq_id",
                "polymer_entity.rcsb_polymer_entity_align.reference_database_name",
                "polymer_entity.rcsb_polymer_entity_align.reference_database_accession",
                "polymer_entity.rcsb_polymer_entity_align.aligned_regions.entity_beg_seq_id",
                "polymer_entity.rcsb_polymer_entity_align.aligned_regions.ref_beg_seq_id",
                "polymer_entity.rcsb_polymer_entity_align.aligned_regions.length",
                "polymer_entity.entry.rcsb_entry_container_identifiers.assembly_ids",
                "rcsb_polymer_instance_info.modeled_residue_count",
                "rcsb_polymer_instance_feature_summary.coverage",
                "rcsb_polymer_instance_feature_summary.type",
                "polymer_entity.rcsb_polymer_entity.pdbx_mutation",
                "polymer_entity.entity_poly.rcsb_mutation_count",
            ]
        )
        r = query.exec()
        out = {}
        for instance in r['data']['polymer_entity_instances']:
            inst_id = instance['rcsb_id']
            bas_id  = instance['polymer_entity']['entry']['rcsb_entry_container_identifiers']['assembly_ids']
            missing_ranges = []
            for feat in instance.get("rcsb_polymer_instance_feature", []):
                if feat.get("type") == "UNOBSERVED_RESIDUE_XYZ":
                    for seg in feat.get("feature_positions", []):
                        missing_ranges.append((seg["beg_seq_id"], seg["end_seq_id"]))
            mapped = {}
            align = instance["polymer_entity"]["rcsb_polymer_entity_align"]
            for entity in align:
                ref_db = entity["reference_database_name"]
                ref_acc = entity["reference_database_accession"]
                if ref_db == 'UniProt' and ref_acc == uniprotid:
                    blocks = entity["aligned_regions"]
                    mapped = {(beg, end): self.pdb2uniprotRange(beg, end, blocks)
                              for (beg, end) in missing_ranges}
                    break
            coverage_field = instance['rcsb_polymer_instance_feature_summary']
            unobs_res_cov = next((d['coverage'] for d in coverage_field if d['type'] == 'UNOBSERVED_RESIDUE_XYZ'), 0)
            modeled_residue_count = instance['rcsb_polymer_instance_info']['modeled_residue_count']
            pdbx_mutation = instance['polymer_entity']['rcsb_polymer_entity']['pdbx_mutation']
            rcsb_mutation_count = instance['polymer_entity']['entity_poly']['rcsb_mutation_count']
            out[inst_id] = {
                'biological_assembly': bas_id,
                'missing_ranges': missing_ranges,
                'mapped': mapped,
                'unobs_res_cov': unobs_res_cov,
                'modeled_residue_count': modeled_residue_count,
                'pdbx_mutation': pdbx_mutation,
                'rcsb_mutation_count': rcsb_mutation_count,
            }
        return out

    def _parseLigandsfromPDB(self, pdblist):
        query = Query(
            input_type="entries",
            input_ids=pdblist,
            return_data_list=[
                "nonpolymer_entities.pdbx_entity_nonpoly.comp_id",
                "nonpolymer_entities.pdbx_entity_nonpoly.name",
                "nonpolymer_entities.rcsb_nonpolymer_entity.pdbx_description",
            ]
        )
        r = query.exec()
        out = {}
        for entry in r['data']['entries']:
            pdbid = entry['rcsb_id']
            ligand_list = entry.get('nonpolymer_entities', None)
            if ligand_list is None:
                out[pdbid] = None
                continue
            out[pdbid] = []
            for entity in ligand_list:
                _dict = {
                    'comp_id': entity['pdbx_entity_nonpoly']['comp_id'],
                    'name': entity['pdbx_entity_nonpoly']['name'],
                    'pdbx_description': entity['rcsb_nonpolymer_entity']['pdbx_description']
                }
                out[pdbid].append(_dict)
        return out

    def _parseRvaluefromPDB(self, pdblist):
        query = Query(
            input_type="entries",
            input_ids=pdblist,
            return_data_list=[
                "rcsb_accession_info.initial_release_date",
                "refine.ls_R_factor_R_free",
                "refine.ls_R_factor_R_work",
                "refine.ls_R_factor_obs",
            ]
        )
        r = query.exec()
        out = {}
        for entry in r['data']['entries']:
            pdbid = entry['rcsb_id']
            rcsb_accession_info = entry.get('rcsb_accession_info', None)
            initial_release_date = rcsb_accession_info['initial_release_date'] if rcsb_accession_info else None
            refine = entry.get('refine', None)
            if refine:
                refine_info = refine[0]
                ls_R_factor_R_free = refine_info['ls_R_factor_R_free']
                ls_R_factor_R_work = refine_info['ls_R_factor_R_work']
                ls_R_factor_obs    = refine_info['ls_R_factor_obs']
            else:
                ls_R_factor_R_free = None
                ls_R_factor_R_work = None
                ls_R_factor_obs    = None
            out[pdbid] = {
                'initial_release_date': initial_release_date,
                'ls_R_factor_R_free': ls_R_factor_R_free,
                'ls_R_factor_R_work': ls_R_factor_R_work,
                'ls_R_factor_obs': ls_R_factor_obs,
            }
        return out

    def _parsePDB(self):
        data = self._rawdata
        PDBdata = {}
        for key, value in data.items():
            if not key.startswith('dbReference'):
                continue
            try:
                pdbid = value['PDB']
            except (KeyError, TypeError):
                continue
            method = value.get('method', None)
            pdbchains = value.get('chains', [])
            resolution = value.get('resolution', '1.00 A')
            resolution = float(resolution.split(' ')[0])
            chains = []
            resrange = None
            try:
                pdbchains_list = comma_splitter(pdbchains)
                for chain in pdbchains_list:
                    chids, resrange = chain.split('=')
                    chids = [chid.strip() for chid in chids.split('/')]
                    for chid in chids:
                        chains.append(chid)
            except Exception as e:
                LOGGER.warn(str(e))
                LOGGER.warn('Suspected no chain information')
            PDBdata[pdbid] = {
                'method': method,
                'resolution': resolution,
                'chains': chains,
                'resrange': resrange,
            }
        pdblist = list(PDBdata.keys())
        if not pdblist:
            self._pdbdata = PDBdata
            return
        ligands = self._parseLigandsfromPDB(pdblist)
        rvalues = self._parseRvaluefromPDB(pdblist)
        auth2label = self.fetchAsymIDs(pdblist)
        for pdbid in PDBdata:
            PDBdata[pdbid]['ligand'] = ligands.get(pdbid)
            rval = rvalues.get(pdbid, {})
            PDBdata[pdbid]['initial_release_date'] = rval.get('initial_release_date')
            PDBdata[pdbid]['ls_R_factor_R_free'] = rval.get('ls_R_factor_R_free')
            PDBdata[pdbid]['ls_R_factor_R_work'] = rval.get('ls_R_factor_R_work')
            PDBdata[pdbid]['ls_R_factor_obs'] = rval.get('ls_R_factor_obs')
            chains = PDBdata[pdbid]['chains']
            asym_record = auth2label.get(pdbid, [])
            if asym_record:
                chains = [next((rec['label_asym_id'] for rec in asym_record if rec['auth_asym_id'] == chid), chid) for chid in chains]
            pdb_instances = [f'{pdbid}.{chid}' for chid in chains]
            seq_annot = self._parseSeqAnnotfromPDB(pdb_instances) if pdb_instances else None
            PDBdata[pdbid]['seq_annot'] = seq_annot
        self._pdbdata = PDBdata

    def fetchAsymIDs(self, pdblist):
        """
        Convert auth_asym_ids to label_asym_id.
        """
        q = Query(
            input_type="entries",
            input_ids=pdblist,
            return_data_list=[
                "rcsb_id",
                "polymer_entities.polymer_entity_instances.rcsb_id",
                "polymer_entities.polymer_entity_instances.rcsb_polymer_entity_instance_container_identifiers.auth_asym_id",
                "polymer_entities.polymer_entity_instances.rcsb_polymer_entity_instance_container_identifiers.asym_id",
                "polymer_entities.polymer_entity_instances.rcsb_polymer_entity_instance_container_identifiers.entity_id",
            ]
        )
        r = q.exec()
        auth2label = {}
        for entry in r['data']['entries']:
            pdbid = entry['rcsb_id']
            asym_record = []
            for entity in entry['polymer_entities']:
                for inst in entity['polymer_entity_instances']:
                    asym_record.append({
                        'rcsb_id': inst['rcsb_id'],
                        'auth_asym_id': inst['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id'],
                        'label_asym_id': inst['rcsb_polymer_entity_instance_container_identifiers']['asym_id'],
                        'entity_id': inst['rcsb_polymer_entity_instance_container_identifiers']['entity_id'],
                    })
            auth2label[pdbid] = asym_record
        return auth2label

    def pdb2uniprotRange(self, beg, end, blocks):
        """Map an entity [beg,end] segment to UniProt ranges."""
        out = []
        for b in blocks:
            e0 = b["entity_beg_seq_id"]
            r0 = b["ref_beg_seq_id"]
            L  = b["length"]
            e1 = e0 + L - 1
            s = max(beg, e0)
            t = min(end, e1)
            if s <= t:
                rs = r0 + (s - e0)
                rt = r0 + (t - e0)
                out.append((rs, rt))
        return out

    def _parse(self):
        LOGGER.info(f'Parse UniProt information of {self.getAccession()}...')
        LOGGER.timeit('_parse')
        self._parseActiveSite()
        self._parseBindingSite()
        self._parseSite()
        self._parseCofactor()
        self._parseDNAbinding()
        self._parseZincfinger()
        self._parseCellLocation()
        self._parsePTMProcessing()
        self._parseFamilyDomains()
        self._parseInteractions()
        self._parseComplexPortals()
        self._parsePDB()
        LOGGER.report(f'Parsing in %.1fs.', '_parse')
        
def searchUniprot(id):
    """Search Uniprot with *id* and return a :class:`UniprotRecord` containing the results. 
    """
    def _queryUniprot(*args, n_attempts=3, dt=1, **kwargs):
        """
        Redefine prody function to check for no internet connection
        """
        attempt = 0
        while attempt < n_attempts:
            try:
                _ = openURL('http://www.uniprot.org/')
                break
            except:
                LOGGER.info(
                    f'Attempt {attempt} to contact www.uniprot.org failed')
                attempt += 1
                time.sleep((attempt+1)*dt)
        else:
            _ = openURL('http://www.uniprot.org/')
        return queryUniprot(*args, **kwargs)

    data = _queryUniprot(id)
    return UniprotRecord(data)

def ChEBI2ligandID(ChEBI):
    # Using ChEBI
    q1 = AttributeQuery(
        attribute="rcsb_chem_comp_related.resource_name",
        operator="exact_match",
        value="ChEBI"  # can also use "ChEMBL", "DrugBank", or "PubChem"
    )
    q2 = AttributeQuery(
        attribute="rcsb_chem_comp_related.resource_accession_code",
        operator="exact_match",
        value=ChEBI,
    )
    q2 = NestedAttributeQuery(q1, q2)
    r = list(q2(return_type="mol_definition"))
    return r[0]

def filter_essence(entry):
    id = entry.get('id')
    # Functional sites
    binding_site = entry.get('binding_site')
    binding_site = [p.get('position') for p in binding_site] if len(binding_site) > 0 else []
    
    active_site = entry.get('active_site')
    active_site = [p.get('position') for p in active_site] if len(active_site) > 0 else []
    
    site = entry.get('site')
    site = [p.get('position') for p in site] if len(site) > 0 else []
    
    dna_binding = entry.get('dna_binding')
    dna_binding = [f"{pos['begin']}-{pos['end']}" for pos in dna_binding] if len(dna_binding) > 0 else []
    
    zinc_finger = entry.get('zinc_finger')
    zinc_finger = [f"{pos['begin']}-{pos['end']}" for pos in zinc_finger] if len(zinc_finger) > 0 else []
    
    # Cellular location
    cell_location = entry['cell_location']
    _intra_mem = []
    _topol_dom = []
    _trans_mem = []
    for elem in cell_location:
        if elem.get('type') == 'intramembrane region':
            _intra_mem.append(f"{elem.get('begin')}-{elem.get('end')}")
        elif elem.get('type') == 'topological domain':
            _topol_dom.append(f"{elem.get('begin')}-{elem.get('end')}")
        elif elem.get('type') == 'transmembrane region':
            _trans_mem.append(f"{elem.get('begin')}-{elem.get('end')}")
        
    return {
        'id': id,
        'protein': entry.get('protein').get('recommend_name'),
        'gene': entry.get('gene'),
        'organism': entry.get('organism').get('common_name'),
        'sequence': entry.get('sequence'),
        'uni_cofactor': entry['cofactor'],
        'intra_mem': _intra_mem,
        'topol_dom': _topol_dom,
        'trans_mem': _trans_mem,
        'binding_site': binding_site,
        'active_site': active_site,
        'site': site,
        'dna_binding': dna_binding,
        'zinc_finger': zinc_finger,
        'pdb': entry.get('pdb'),
        'alphafold': entry.get('alphafold'),
    }

def pickPDBfromUniprot(entry, to_file=None):
    data = filter_essence(entry)
    id = data['id'] # uniprot id
    seq_len = len(data['sequence']) # uniprot protein length
    
    # Cofactor
    uni_cofactor = data['uni_cofactor']
    uni_cofactor_chebi = [cf['chebi'] for cf in uni_cofactor]
    uni_cofactor_id = []
    for l in uni_cofactor_chebi:
        uni_cofactor_id.extend(COFACTOR[l])
    if len(uni_cofactor_id) == 0:
        LOGGER.info(f"{id}, no cofactor")
    else: 
        LOGGER.info(f"{id}, cofactor(s):")
        for i, cof in enumerate(uni_cofactor):
            LOGGER.info(f"{cof['name']}, {cof['chebi']}, {uni_cofactor_id[i]}")
    
    # Functional site
    _func_site_ = []
    for sname in ['active_site', 'binding_site', 'dna_binding', 'zinc_finger']:
        _site = data[sname]
        if len(_site) != 0:
            _func_site_.extend(_site)
            LOGGER.info(f"{sname}: {_site}")
    
    # lst = ['280', '99-101', '103-105']
    # out: [280, 99, 100, 101, 103, 104, 105]
    func_site = []
    for item in _func_site_:
        item = str(item)
        if '-' in item:
            start, end = map(int, item.split('-'))
            func_site.extend(range(start, end + 1))
        else:
            func_site.append(int(item))
    
    # PDB structures
    pdb = data.get('pdb')
    if len(pdb) == 0:
        LOGGER.warn(f"{id} has no PDB structures")
    
    # Retrieve key information to rank PDB structures
    PDBrank = []
    for pdbid, entry in pdb.items():
        resolution = entry['resolution']
        pdb_ligand = entry['ligand']
        r_free = entry['ls_R_factor_R_free']
        seq_annot = entry['seq_annot']
        resrange = entry['resrange']
        
        # Retrieve ligand id
        if pdb_ligand is not None:
            pdb_ligand_id = [l['comp_id'] for l in pdb_ligand]
        else:
            pdb_ligand_id = []
            
        # Count number of _cofactors = cofactor coenzyme + cofactor ion
        # Count number of (ligands+drugs) other than _cofactors
        _cofactors = set(pdb_ligand_id).intersection(set(uni_cofactor_id))
        _coenzymes = set(pdb_ligand_id).intersection(set(COFACTOR_COENZYME))
        _cofactors = _cofactors.union(_coenzymes)
        _ligands_or_drugs = set(pdb_ligand_id).difference(set(uni_cofactor_id))
        
        # resolved_len with respect to uniprot sequence e.g., '100-250'
        if resrange is None:
            resolved_len = 0
            LOGGER.warn(f"{pdbid} has no resolved range")
        else:
            _split = resrange.split('-')
            resolved_len = int(_split[1]) - int(_split[0]) + 1
        
        # Iterate each chain (instance) recorded in uniprot
        for inst_id, value in seq_annot.items():
            n_mut = value['rcsb_mutation_count']
            bas_id = value['biological_assembly']
            
            # dict, key: missing pdb resID, value: uniprot ID 
            # e.g., (1, 12): [(24, 30)] 
            mapped = value['mapped'] 
            
            # modeled_len is resolved_len substracted missing residues
            n_missing_res = 0
            n_missing_fsite = 0
            for _, m in mapped.items():
                if len(m) == 0:
                    continue
                _range = m[0]
                missing_range = range(_range[0], _range[1])
                missing_fsite = set(func_site).intersection(set(missing_range))
                n_missing_res += len(missing_range)
                n_missing_fsite += len(missing_fsite)
                
            modeled_len = resolved_len - n_missing_res
            coverage = modeled_len / seq_len
            PDBrank.append(
                (
                    id, uni_cofactor_id, seq_len, # 0, 1, 2
                    pdbid, inst_id, bas_id, pdb_ligand_id, # 3, 4, 5, 6
                    _cofactors, len(_cofactors), # 7, 8
                    _ligands_or_drugs, len(_ligands_or_drugs), # 9, 10
                    resrange, modeled_len, n_mut, n_missing_fsite, coverage, # 11, 12, 13, 14, 15
                    resolution, r_free # 16, 17
                )
            )
    # len(_cofactors) → n_missing_fsite → n_mut →
    # coverage → resolution → r_free →
    # len(_ligands_or_drugs)
    PDBrank.sort(key=lambda x: (-x[8], -x[14], x[13], 
                                -x[15], x[16], x[17], 
                                x[10])) # Smallest score first    
    # Save to csv file
    if to_file is not None:
        import pandas as pd 
        columns = [
            'id', 'uni_cofactor', 'uni_seq_len',
            'pdbid', 'chid', 'basid', 'all_pdb_ligands',
            'pdb_cofactors', 'n_pdb_cofactors', 
            'pdb_ligand_or_drugs', 'n_pdb_ligand_or_drugs',
            'resrange', 'modeled_len', 'n_mut', 'n_missing_fsite', 'coverage', 
            'resolution', 'r_free'
        ]
        df = pd.DataFrame(columns=columns, data=PDBrank)
        df.to_csv(to_file, index=False)
    return PDBrank


from rcsbapi.search import AttributeQuery, NestedAttributeQuery
from rcsbapi.sequence import Alignments
from src.utils.logger import LOGGER 

def queryUid2RCSB(uid):
    """Query UniProtID to collect all PDB ID (entities and instances)"""
    # uniprot_id = 'P0DMS8'
    q_acc = AttributeQuery(
        attribute='rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession',
        value=uid,
        operator="exact_match",
    )
    q_db = AttributeQuery(
        attribute='rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name',
        value="UniProt",
        operator="exact_match",
    )
    # Group nested attributes using `NestedAttributeQuery`
    query = NestedAttributeQuery(q_acc, q_db)
    
    # Query data
    pdb_id = list(query())
    pdb_entity = list(query(return_type='polymer_entity'))
    pdb_instance = list(query(return_type='polymer_instance'))

    # Retrieve
    out = []
    for id in pdb_id:
        _entity = next(
            (_ent for _ent in pdb_entity if _ent.startswith(f"{id}_")), None
        )
        if not _entity:
            continue
        
        _instance = [_ins for _ins in pdb_instance if _ins.startswith(f"{id}.")]
        sample = {}
        sample["pdb_id"] = _entity[:4]
        sample["pdb_entity"] = _entity
        sample["pdb_instance"] = _instance
        out.append(sample)
    return out

def align_RCSB2UniProt(pdb_entity, uid):

    # Fetch alignments between a UniProt Accession and PDB Entity 
    # (i.e., unique sequence within the entry)
    query = Alignments(
        db_from="PDB_ENTITY",
        db_to="UNIPROT",
        query_id=pdb_entity,
        return_data_list=["query_sequence", "target_alignments"]
    )
    r = query.exec()
    
    uid_record = next(
        (item for item in r['data']['alignments']['target_alignments'] \
        if item['target_id'] == uid), None
    )
    if not uid_record:
        LOGGER.warn(f"No information for {pdb_entity}")
        return None

    aligned_regions = uid_record['aligned_regions']
    query_seq = r['data']['alignments']['query_sequence']
    target_seq = uid_record['target_sequence']

    if len(aligned_regions) > 1:
        LOGGER.warn(f"More than 1 aligned regions for {pdb_entity} and {uid}")

    for i, region in enumerate(aligned_regions):
        pdb_seq = query_seq[region['query_begin']-1 : region['query_end']]
        uniprot_seq = target_seq[region['target_begin']-1 : region['target_end']]

        # Find residues that disagree between PDB and UniProt 
        difference = []
        for j, (pdb_aa, uniprot_aa) in enumerate(zip(pdb_seq, uniprot_seq)):
            if pdb_aa != uniprot_aa:
                difference.append({
                    'pdb_resID': region['query_begin']+j,
                    'uniprot_resID': region['target_begin']+j,
                    'uniprot_aa': uniprot_aa,
                    'pdb_aa': pdb_aa,
                })
        aligned_regions[i]['difference'] = difference
    return {
        'aligned_regions': aligned_regions,
        'query_seq': query_seq,
        'target_seq': target_seq,   
    }

# uid = 'P02925'
# data = queryUid2RCSB(uid)
# for data_idx in range(len(data)):
#     sample = data[data_idx]
#     pdb_entity = sample['pdb_entity']

#     align = align_RCSB2UniProt(pdb_entity, uid)
#     if not align:
#         continue

#     sample['aligned_regions'] = align['aligned_regions']
#     sample['query_seq'] = align['query_seq']
#     sample['target_seq'] = align['target_seq']