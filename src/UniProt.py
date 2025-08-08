# -*- coding: utf-8 -*-
"""This module defines a class and relative functions for mapping Uniprot
sequences to PDB and Pfam databases."""

import re
import time
import re

from prody import queryUniprot
from prody.utilities import openURL
from rcsbapi.data import DataQuery as Query

from .utils.logger import LOGGER

comma_splitter = re.compile(r'\s*,\s*').split
ns = {'up': 'http://uniprot.org/uniprot'}
    
# -*- coding: utf-8 -*-
"""This module defines a class and relative functions for mapping Uniprot
sequences to PDB and Pfam databases."""
from rcsbapi.data import DataQuery as Query
import os
import re
import dill as pickle
import datetime
import time
import numpy as np
import urllib.parse
import requests 
import re
import traceback

import prody
from prody import parsePDB, Atomic, queryUniprot
from prody.utilities import openURL
from Bio.pairwise2 import align as bioalign
from Bio.pairwise2 import format_alignment

from prody import LOGGER

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

comma_splitter = re.compile(r'\s*,\s*').split
ns = {'up': 'http://uniprot.org/uniprot'}
    
class UniprotRecord(object):
    """This class provides a wrapper for UniProt data including functions 
    for accessing particular fields and parsing associated PDB entries."""
    def __init__(self, data):
        self._rawdata = data
        self._pdbdata = []
        self._parse()

    def __repr__(self):
        return '<UniprotRecord: %s>'%self.getTitle()

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
        """accession tag"""
        return self.getEntry('accession', index)
    
    def getName(self, index=0):
        """name tag"""
        return self.getEntry('name', index)

    def getProtein(self, index=0):
        """protein tag
        <protein>
            <recommendedName>
                <fullName>Gap junction beta-2 protein</fullName>
            </recommendedName>
            <alternativeName>
                <fullName evidence="57 58">Connexin-26</fullName>
                <shortName evidence="58">Cx26</shortName>
            </alternativeName>
        </protein>
        """
        protein = self.getEntry('protein', index)
        
        try:
            recommend_elem = protein.find('up:recommendedName', ns)
            alternative_elem = protein.find('up:recommendedName', ns)
            
            recommend_name = recommend_elem.find('up:fullName', ns)
            alter_fullname = alternative_elem.find('up:fullName', ns)
            alter_shortname = alternative_elem.find('up:shortName', ns)
            
            recommend_name = recommend_name.text if recommend_name is not None else None
            alter_fullname = alter_fullname.text if alter_fullname is not None else None
            alter_shortname = alter_shortname.text if alter_shortname is not None else None
        except:
            submitted_name = protein.find('up:submittedName/up:fullName', ns)
            submitted_name = submitted_name.text if submitted_name is not None else None
            recommend_name = submitted_name
            alter_fullname = None
            alter_shortname = None
        
        return {
            'recommend_name': recommend_name,
            'alter_fullname': alter_fullname,
            'alter_shortname': alter_shortname            
        }
        
    def getGene(self, index=0):
        """gene tag
        <gene>
            <name type="primary">GJB2</name>
        </gene>
        """
        try:
            gene = self.getEntry('gene', index)
            name_elem = gene.find('up:name[@type="primary"]', ns)
            return name_elem.text if name_elem is not None else None
        
        except Exception as e:
            LOGGER.warn(f'Error while parsing {id}: {e} -> None')
            return None

    def getOrganism(self, index=0):
        """organism tag
        <organism>
            <name type="scientific">Homo sapiens</name>
            <name type="common">Human</name>
            <dbReference type="NCBI Taxonomy" id="9606"/>
            <lineage>
                <taxon>Eukaryota</taxon>
                <taxon>Metazoa</taxon>
                <taxon>Chordata</taxon>
                <taxon>Craniata</taxon>
                <taxon>Vertebrata</taxon>
                <taxon>Euteleostomi</taxon>
                <taxon>Mammalia</taxon>
                <taxon>Eutheria</taxon>
                <taxon>Euarchontoglires</taxon>
                <taxon>Primates</taxon>
                <taxon>Haplorrhini</taxon>
                <taxon>Catarrhini</taxon>
                <taxon>Hominidae</taxon>
                <taxon>Homo</taxon>
            </lineage>
        </organism>
        """
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
        
    def getReference(self, index=0):
        """reference tag"""
        pass
    
    def getComment(self, index=0):
        """conmment tag"""
        pass
        
    def getDBreference(self, index=0):
        """dbReference tag"""
        pass
    
    def getProteinExistence(self, index=0):
        """proteinExistence tag"""
        pass
     
    def getKeyword(self, index=0):
        """keyword tag"""
        pass
    
    def getFeature(self, index=0):
        """feature tag"""
        pass
    
    def getEvidence(self, index=0):
        """evidence tag"""
        pass
    
    def getSequence(self, index=0):
        return self.getEntry('sequence', index)
    
    def getZincFinger(self):
        return self._zinc_finger
    
    def getDNAbinding(self):
        return self._dna_binding
    
    def getActivateSite(self):
        return self._active_site
    
    def getBindingSite(self):
        return self._binding_site
    
    def getSite(self):
        return self._site
    
    def getAlphaFold(self):
        """<dbReference type="AlphaFoldDB" id="Q9UDY8"/>"""
        AlphaFoldDB = None
        for key, value in self._rawdata.items():
            if not key.startswith('dbReference'):
                continue

            if type(value) != list or len(value) != 2:
                continue
            
            # [('type', 'AlphaFoldDB'), ('id', 'Q13509')]
            if value[0][1] == 'AlphaFoldDB':
                AlphaFoldDB = value[1][1]
                break
        return AlphaFoldDB
    
    def getCofactor(self):
        return self._cofactors
    
    def getTitle(self):
        uid = self.getAccession()
        name = self.getName()
        return '%s (%s)'%(uid, name)

    def getEntry(self, item, index=0):
        key = '%s%4d'%(item, index)
        if key in self._rawdata:
            return self._rawdata[key]
        else:
            raise KeyError('%s does not exist in the Uniprot record'%key)

    def _parseDNAbinding(self):
        data = self._rawdata
        dna_binding = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            if value.get('type') != "DNA-binding region":
                continue
            """
            <feature type="DNA-binding region" description="HMG box 1" evidence="4">
                <location>
                <begin position="9"/>
                <end position="79"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            dna_binding.append({
                'description': descp, 
                'begin': begin, 
                'end': end
            })
        self._dna_binding = dna_binding
            
    def _parseZincfinger(self):
        data = self._rawdata
        zinc_finger = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            if value.get('type') != "zinc finger region":
                continue
            """
            <feature type="zinc finger region" description="C2H2-type 1" evidence="1">
                <location>
                <begin position="110"/>
                <end position="133"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            zinc_finger.append({
                'description': descp, 
                'begin': begin, 
                'end': end
            })
        self._zinc_finger = zinc_finger

    def _parseActiveSite(self):
        data = self._rawdata
        active_site = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            
            if value.get('type') != "active site":
                continue
            """
            <feature type="active site" description="Proton donor" evidence="2">
                <location>
                <position position="613"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            pos_elem = value.find('up:location/up:position', ns)
            pos   = int(pos_elem.attrib.get('position')) if pos_elem is not None else None
            active_site.append({
                'description': descp, 
                'position': pos
            })
        self._active_site = active_site
    
    def _parseBindingSite(self):
        data = self._rawdata
        binding_site = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            
            if value.get('type') != "binding site":
                continue
            
            """
            <feature type="binding site" evidence="7 9 22 23 24">
                <location>
                <position position="617"/>
                </location>
                <ligand>
                <name>Zn(2+)</name>
                <dbReference type="ChEBI" id="CHEBI:29105"/>
                </ligand>
            </feature>
            """
            pos_elem = value.find('up:location/up:position', ns)
            pos   = int(pos_elem.attrib.get('position')) if pos_elem is not None else None
            
            ligand_elem = value.find('up:ligand', ns)
            ligand_name = ligand_elem.find('up:name', ns)
            ligand_name = ligand_name.text if ligand_name is not None else None
            ligand_chebi= ligand_elem.find('up:dbReference[@type="ChEBI"]', ns)
            ligand_chebi = ligand_chebi.attrib['id'] if ligand_chebi is not None else None
            binding_site.append({
                'position': pos, 
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
            
            """
            <feature type="site" description="Breakpoint for translocation to form BIRC2-MALT1">
                <location>
                    <begin position="323"/>
                    <end position="324"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            site.append({
                'description': descp, 
                'begin': begin, 
                'end': end
            })
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
            """
            <comment type="cofactor">
                <cofactor evidence="2">
                    <name>cf_name</name>
                    <dbReference type="ChEBI" id=cf_chebi/>
                </cofactor>
            </comment>
            """
            cf_elem = value.find('up:cofactor', ns)
            # ---
            cf_name = cf_elem.find('up:name', ns)
            cf_chebi= cf_elem.find('up:dbReference[@type="ChEBI"]', ns)
            cf_chebi = cf_chebi.attrib['id'] if cf_chebi is not None else None
            # ---
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
            
            type_list = ['topological domain', 'transmembrane region', 'intramembrane region']
            loc_type  = value.get('type')
            if loc_type not in type_list:
                continue
            """
            <feature type="intramembrane region" evidence="45">
                <location>
                <begin position="2"/>
                <end position="13"/>
                </location>
            </feature>
            """
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

    def _parseSeqAnnotfromPDB(self, pdb_instances):
        """pdb_instances ['1AID.A', '1AID.B']
        """
        query = Query(
            input_type="polymer_entity_instances",
            input_ids=pdb_instances,
            return_data_list=[
                "polymer_entity_instances.rcsb_id",
                "rcsb_polymer_instance_feature.type",
                "rcsb_polymer_instance_feature.feature_positions.beg_seq_id",
                "rcsb_polymer_instance_feature.feature_positions.end_seq_id",
                
                'rcsb_polymer_instance_info.modeled_residue_count',
                'rcsb_polymer_instance_feature_summary.coverage',
                'rcsb_polymer_instance_feature_summary.type',
                'polymer_entity.rcsb_polymer_entity.pdbx_mutation',
                'polymer_entity.entity_poly.rcsb_mutation_count',
            ]
        )
        r = query.exec()

        out = {}
        for entry in r['data']['polymer_entity_instances']:
            instance = entry['rcsb_id']
            instance_feature = entry.get('rcsb_polymer_instance_feature', None)
                
            unobs_res = [ele['feature_positions'] for ele in instance_feature if ele['type'] == 'UNOBSERVED_RESIDUE_XYZ']
            unobs_atom= [ele['feature_positions'] for ele in instance_feature if ele['type'] == 'UNOBSERVED_ATOM_XYZ']
            
            # Turn into a list of residues / atoms
            unobs_res = [
                i for block in unobs_res for r in block
                for i in range(r['beg_seq_id'], r['end_seq_id'] + 1)
            ]
            unobs_atom = [
                i for block in unobs_atom for r in block
                for i in range(r['beg_seq_id'], r['end_seq_id'] + 1)
            ]
            
            # Extract specific types
            coverage_field = entry['rcsb_polymer_instance_feature_summary']    
            unobs_res_cov = next((d['coverage'] for d in coverage_field if d['type'] == 'UNOBSERVED_RESIDUE_XYZ'), 0)
            unobs_atom_cov = next((d['coverage'] for d in coverage_field if d['type'] == 'UNOBSERVED_ATOM_XYZ'), 0)

            coverage = {
                "unobs_res": unobs_res_cov,
                "unobs_atom": unobs_atom_cov
            }
            
            modeled_residue_count = entry['rcsb_polymer_instance_info']['modeled_residue_count']
            pdbx_mutation = entry['polymer_entity']['rcsb_polymer_entity']['pdbx_mutation']
            rcsb_mutation_count = entry['polymer_entity']['entity_poly']['rcsb_mutation_count']
            
            out[instance] = {
                'unobs_res': unobs_res,
                'unobs_atom': unobs_atom,
                'coverage': coverage,
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
            # Return None if no ligand
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
            # Return None if no ligand
            if rcsb_accession_info is None:
                initial_release_date = None
            else:
                initial_release_date = rcsb_accession_info['initial_release_date']
                
            refine = entry.get('refine', None)
            # Return None if no refine
            if refine is None:
                ls_R_factor_R_free = None
                ls_R_factor_R_work = None
                ls_R_factor_obs = None
            else:
                refine_info = refine[0]
                ls_R_factor_R_free = refine_info['ls_R_factor_R_free']
                ls_R_factor_R_work = refine_info['ls_R_factor_R_work']
                ls_R_factor_obs    = refine_info['ls_R_factor_obs']
                
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
            except (KeyError, TypeError) as e:
                continue
            """
            <dbReference type="PDB" id=pdbid>
                <property type="method" value="EM"/>
                <property type="resolution" value=resolution/>
                <property type="chains" value=pdbchains/>
    		</dbReference>
            """
            method = value['method']
            method = value.get('method', None)
            # pdbchains = value['chains'] # e.g. "B/D/F/G/H/I=1-450"
            pdbchains = value.get('chains', []) # e.g. "B/D/F/G/H/I=1-450"
            resolution = value.get('resolution', '1.00 A')
            resolution = float(resolution.split(' ')[0])
            
            # example chain strings: "A=27-139, B=140-150" or "A/B=27-150"
            chains = []
            resrange = None
            try:
                pdbchains = comma_splitter(pdbchains)
                for chain in pdbchains:
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
        if len(pdblist) == 0:
            self._pdbdata = PDBdata
            return
        
        # RCSB Data API: entries 
        # Retrieved info: ligands, released date, Observed Residual factor (R-value obs) 
        ligands = self._parseLigandsfromPDB(pdblist)
        rvalues = self._parseRvaluefromPDB(pdblist)
        
        # fetchAsymIDs to convert auth_asym_ids into label_asym_id
        auth2label = self.fetchAsymIDs(pdblist)
        
        for pdbid in PDBdata:
            PDBdata[pdbid]['ligand'] = ligands[pdbid]
            PDBdata[pdbid]['initial_release_date'] = rvalues[pdbid]['initial_release_date']
            PDBdata[pdbid]['ls_R_factor_R_free'] = rvalues[pdbid]['ls_R_factor_R_free']
            PDBdata[pdbid]['ls_R_factor_R_work'] = rvalues[pdbid]['ls_R_factor_R_work']
            PDBdata[pdbid]['ls_R_factor_obs'] = rvalues[pdbid]['ls_R_factor_obs']
            
            # fetchAsymIDs to convert auth_asym_ids into label_asym_id
            chains = PDBdata[pdbid]['chains']
            chain_dict = auth2label.get(pdbid, None)
            if chain_dict is None:
                LOGGER.warn(f'fetchAsymIDs: No infor. of {pdbid}')
            else:
                chains = [chain_dict.get(chid, chid) for chid in chains]
            
            # RCSB Data API: polymer_entity_instances 
            # Retrieved info: Sequence Annotations - UNOBSERVED_RESIDUE_XYZ
            pdb_instances = [f'{pdbid}.{chid}' for chid in chains]
            if len(pdb_instances) > 0:
                seq_annot = self._parseSeqAnnotfromPDB(pdb_instances)
            else:
                seq_annot = None
            PDBdata[pdbid]['seq_annot'] = seq_annot
        self._pdbdata = PDBdata
    
    def fetchAsymIDs(self, pdblist):
        """
        Convert auth_asym_ids into label_asym_id. 
        https://www.rcsb.org/docs/general-help/identifiers-in-pdb#:~:text=type%20of%20entity.-,Macromolecular%20Instance%20ID,R%2C%20while%20the%20PDB%20assigned%20ones%20are%20C%20and%20D%20respectively.,-The%20polymer%20sequences
        
        Return: dict
            given auth_asym_ids, dict can be used as a look up table to find label_asym_id.
        """
        query = Query(
            input_type="entries",
            input_ids=pdblist,
            return_data_list=[
                "polymer_entities.polymer_entity_instances.rcsb_id",
                "polymer_entities.rcsb_polymer_entity_container_identifiers.auth_asym_ids",
            ]
        )
        r = query.exec()

        if len(r['data']['entries']) == 0:
            LOGGER.warn(f'fetchAsymIDs: Check input {pdblist}')
            return {}

        auth2label = {}
        for entry in r['data']['entries']:
            pdbid = entry['rcsb_id']
            
            _dict = {}
            for entity in entry['polymer_entities']:
                label_asym_id = entity['polymer_entity_instances'][0]['rcsb_id'].split('.')[-1]
                auth_asym_id  = entity['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'][0]
                _dict[auth_asym_id] = label_asym_id
            
            auth2label[pdbid] = _dict
        return auth2label

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

comma_splitter = re.compile(r'\s*,\s*').split
ns = {'up': 'http://uniprot.org/uniprot'}
    
class UniprotRecord(object):
    """This class provides a wrapper for UniProt data including functions 
    for accessing particular fields and parsing associated PDB entries."""
    def __init__(self, data):
        self._rawdata = data
        self._pdbdata = []
        self._parse()

    def __repr__(self):
        return '<UniprotRecord: %s>'%self.getTitle()

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
        """accession tag"""
        return self.getEntry('accession', index)
    
    def getName(self, index=0):
        """name tag"""
        return self.getEntry('name', index)

    def getProtein(self, index=0):
        """protein tag
        <protein>
            <recommendedName>
                <fullName>Gap junction beta-2 protein</fullName>
            </recommendedName>
            <alternativeName>
                <fullName evidence="57 58">Connexin-26</fullName>
                <shortName evidence="58">Cx26</shortName>
            </alternativeName>
        </protein>
        """
        protein = self.getEntry('protein', index)
        
        try:
            recommend_elem = protein.find('up:recommendedName', ns)
            alternative_elem = protein.find('up:recommendedName', ns)
            
            recommend_name = recommend_elem.find('up:fullName', ns)
            alter_fullname = alternative_elem.find('up:fullName', ns)
            alter_shortname = alternative_elem.find('up:shortName', ns)
            
            recommend_name = recommend_name.text if recommend_name is not None else None
            alter_fullname = alter_fullname.text if alter_fullname is not None else None
            alter_shortname = alter_shortname.text if alter_shortname is not None else None
        except:
            submitted_name = protein.find('up:submittedName/up:fullName', ns)
            submitted_name = submitted_name.text if submitted_name is not None else None
            recommend_name = submitted_name
            alter_fullname = None
            alter_shortname = None
        
        return {
            'recommend_name': recommend_name,
            'alter_fullname': alter_fullname,
            'alter_shortname': alter_shortname            
        }
        
    def getGene(self, index=0):
        """gene tag
        <gene>
            <name type="primary">GJB2</name>
        </gene>
        """
        try:
            gene = self.getEntry('gene', index)
            name_elem = gene.find('up:name[@type="primary"]', ns)
            return name_elem.text if name_elem is not None else None
        
        except Exception as e:
            LOGGER.warn(f'Error while parsing {id}: {e} -> None')
            return None

    def getOrganism(self, index=0):
        """organism tag
        <organism>
            <name type="scientific">Homo sapiens</name>
            <name type="common">Human</name>
            <dbReference type="NCBI Taxonomy" id="9606"/>
            <lineage>
                <taxon>Eukaryota</taxon>
                <taxon>Metazoa</taxon>
                <taxon>Chordata</taxon>
                <taxon>Craniata</taxon>
                <taxon>Vertebrata</taxon>
                <taxon>Euteleostomi</taxon>
                <taxon>Mammalia</taxon>
                <taxon>Eutheria</taxon>
                <taxon>Euarchontoglires</taxon>
                <taxon>Primates</taxon>
                <taxon>Haplorrhini</taxon>
                <taxon>Catarrhini</taxon>
                <taxon>Hominidae</taxon>
                <taxon>Homo</taxon>
            </lineage>
        </organism>
        """
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
        
    def getReference(self, index=0):
        """reference tag"""
        pass
    
    def getComment(self, index=0):
        """conmment tag"""
        pass
        
    def getDBreference(self, index=0):
        """dbReference tag"""
        pass
    
    def getProteinExistence(self, index=0):
        """proteinExistence tag"""
        pass
     
    def getKeyword(self, index=0):
        """keyword tag"""
        pass
    
    def getFeature(self, index=0):
        """feature tag"""
        pass
    
    def getEvidence(self, index=0):
        """evidence tag"""
        pass
    
    def getSequence(self, index=0):
        return self.getEntry('sequence', index)
    
    def getZincFinger(self):
        return self._zinc_finger
    
    def getDNAbinding(self):
        return self._dna_binding
    
    def getActivateSite(self):
        return self._active_site
    
    def getBindingSite(self):
        return self._binding_site
    
    def getSite(self):
        return self._site
    
    def getAlphaFold(self):
        """<dbReference type="AlphaFoldDB" id="Q9UDY8"/>"""
        AlphaFoldDB = None
        for key, value in self._rawdata.items():
            if not key.startswith('dbReference'):
                continue

            if type(value) != list or len(value) != 2:
                continue
            
            # [('type', 'AlphaFoldDB'), ('id', 'Q13509')]
            if value[0][1] == 'AlphaFoldDB':
                AlphaFoldDB = value[1][1]
                break
        return AlphaFoldDB
    
    def getCofactor(self):
        return self._cofactors
    
    def getTitle(self):
        uid = self.getAccession()
        name = self.getName()
        return '%s (%s)'%(uid, name)

    def getEntry(self, item, index=0):
        key = '%s%4d'%(item, index)
        if key in self._rawdata:
            return self._rawdata[key]
        else:
            raise KeyError('%s does not exist in the Uniprot record'%key)

    def _parseDNAbinding(self):
        data = self._rawdata
        dna_binding = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            if value.get('type') != "DNA-binding region":
                continue
            """
            <feature type="DNA-binding region" description="HMG box 1" evidence="4">
                <location>
                <begin position="9"/>
                <end position="79"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            dna_binding.append({
                'description': descp, 
                'begin': begin, 
                'end': end
            })
        self._dna_binding = dna_binding
            
    def _parseZincfinger(self):
        data = self._rawdata
        zinc_finger = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            if value.get('type') != "zinc finger region":
                continue
            """
            <feature type="zinc finger region" description="C2H2-type 1" evidence="1">
                <location>
                <begin position="110"/>
                <end position="133"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            zinc_finger.append({
                'description': descp, 
                'begin': begin, 
                'end': end
            })
        self._zinc_finger = zinc_finger

    def _parseActiveSite(self):
        data = self._rawdata
        active_site = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            
            if value.get('type') != "active site":
                continue
            """
            <feature type="active site" description="Proton donor" evidence="2">
                <location>
                <position position="613"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            pos_elem = value.find('up:location/up:position', ns)
            pos   = int(pos_elem.attrib.get('position')) if pos_elem is not None else None
            active_site.append({
                'description': descp, 
                'position': pos
            })
        self._active_site = active_site
    
    def _parseBindingSite(self):
        data = self._rawdata
        binding_site = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            
            if value.get('type') != "binding site":
                continue
            
            """
            <feature type="binding site" evidence="7 9 22 23 24">
                <location>
                <position position="617"/>
                </location>
                <ligand>
                <name>Zn(2+)</name>
                <dbReference type="ChEBI" id="CHEBI:29105"/>
                </ligand>
            </feature>
            """
            pos_elem = value.find('up:location/up:position', ns)
            pos   = int(pos_elem.attrib.get('position')) if pos_elem is not None else None
            
            ligand_elem = value.find('up:ligand', ns)
            ligand_name = ligand_elem.find('up:name', ns)
            ligand_name = ligand_name.text if ligand_name is not None else None
            ligand_chebi= ligand_elem.find('up:dbReference[@type="ChEBI"]', ns)
            ligand_chebi = ligand_chebi.attrib['id'] if ligand_chebi is not None else None
            binding_site.append({
                'position': pos, 
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
            
            """
            <feature type="site" description="Breakpoint for translocation to form BIRC2-MALT1">
                <location>
                    <begin position="323"/>
                    <end position="324"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            site.append({
                'description': descp, 
                'begin': begin, 
                'end': end
            })
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
            """
            <comment type="cofactor">
                <cofactor evidence="2">
                    <name>cf_name</name>
                    <dbReference type="ChEBI" id=cf_chebi/>
                </cofactor>
            </comment>
            """
            cf_elem = value.find('up:cofactor', ns)
            # ---
            cf_name = cf_elem.find('up:name', ns)
            cf_chebi= cf_elem.find('up:dbReference[@type="ChEBI"]', ns)
            cf_chebi = cf_chebi.attrib['id'] if cf_chebi is not None else None
            # ---
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
            
            type_list = ['topological domain', 'transmembrane region', 'intramembrane region']
            loc_type  = value.get('type')
            if loc_type not in type_list:
                continue
            """
            <feature type="intramembrane region" evidence="45">
                <location>
                <begin position="2"/>
                <end position="13"/>
                </location>
            </feature>
            """
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

    def _parseSeqAnnotfromPDB(self, pdb_instances):
        """pdb_instances ['1AID.A', '1AID.B']
        """
        query = Query(
            input_type="polymer_entity_instances",
            input_ids=pdb_instances,
            return_data_list=[
                "polymer_entity_instances.rcsb_id",
                "rcsb_polymer_instance_feature.type",
                "rcsb_polymer_instance_feature.feature_positions.beg_seq_id",
                "rcsb_polymer_instance_feature.feature_positions.end_seq_id",
                
                'rcsb_polymer_instance_info.modeled_residue_count',
                'rcsb_polymer_instance_feature_summary.coverage',
                'rcsb_polymer_instance_feature_summary.type',
                'polymer_entity.rcsb_polymer_entity.pdbx_mutation',
                'polymer_entity.entity_poly.rcsb_mutation_count',
            ]
        )
        r = query.exec()

        out = {}
        for entry in r['data']['polymer_entity_instances']:
            instance = entry['rcsb_id']
            instance_feature = entry.get('rcsb_polymer_instance_feature', None)
                
            unobs_res = [ele['feature_positions'] for ele in instance_feature if ele['type'] == 'UNOBSERVED_RESIDUE_XYZ']
            unobs_atom= [ele['feature_positions'] for ele in instance_feature if ele['type'] == 'UNOBSERVED_ATOM_XYZ']
            
            # Turn into a list of residues / atoms
            unobs_res = [
                i for block in unobs_res for r in block
                for i in range(r['beg_seq_id'], r['end_seq_id'] + 1)
            ]
            unobs_atom = [
                i for block in unobs_atom for r in block
                for i in range(r['beg_seq_id'], r['end_seq_id'] + 1)
            ]
            
            # Extract specific types
            coverage_field = entry['rcsb_polymer_instance_feature_summary']    
            unobs_res_cov = next((d['coverage'] for d in coverage_field if d['type'] == 'UNOBSERVED_RESIDUE_XYZ'), 0)
            unobs_atom_cov = next((d['coverage'] for d in coverage_field if d['type'] == 'UNOBSERVED_ATOM_XYZ'), 0)

            coverage = {
                "unobs_res": unobs_res_cov,
                "unobs_atom": unobs_atom_cov
            }
            
            modeled_residue_count = entry['rcsb_polymer_instance_info']['modeled_residue_count']
            pdbx_mutation = entry['polymer_entity']['rcsb_polymer_entity']['pdbx_mutation']
            rcsb_mutation_count = entry['polymer_entity']['entity_poly']['rcsb_mutation_count']
            
            out[instance] = {
                'unobs_res': unobs_res,
                'unobs_atom': unobs_atom,
                'coverage': coverage,
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
            # Return None if no ligand
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
            # Return None if no ligand
            if rcsb_accession_info is None:
                initial_release_date = None
            else:
                initial_release_date = rcsb_accession_info['initial_release_date']
                
            refine = entry.get('refine', None)
            # Return None if no refine
            if refine is None:
                ls_R_factor_R_free = None
                ls_R_factor_R_work = None
                ls_R_factor_obs = None
            else:
                refine_info = refine[0]
                ls_R_factor_R_free = refine_info['ls_R_factor_R_free']
                ls_R_factor_R_work = refine_info['ls_R_factor_R_work']
                ls_R_factor_obs    = refine_info['ls_R_factor_obs']
                
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
            except (KeyError, TypeError) as e:
                continue
            """
            <dbReference type="PDB" id=pdbid>
                <property type="method" value="EM"/>
                <property type="resolution" value=resolution/>
                <property type="chains" value=pdbchains/>
    		</dbReference>
            """
            method = value['method']
            method = value.get('method', None)
            # pdbchains = value['chains'] # e.g. "B/D/F/G/H/I=1-450"
            pdbchains = value.get('chains', []) # e.g. "B/D/F/G/H/I=1-450"
            resolution = value.get('resolution', '1.00 A')
            resolution = float(resolution.split(' ')[0])
            
            # example chain strings: "A=27-139, B=140-150" or "A/B=27-150"
            chains = []
            resrange = None
            try:
                pdbchains = comma_splitter(pdbchains)
                for chain in pdbchains:
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
        if len(pdblist) == 0:
            self._pdbdata = PDBdata
            return
        
        # RCSB Data API: entries 
        # Retrieved info: ligands, released date, Observed Residual factor (R-value obs) 
        ligands = self._parseLigandsfromPDB(pdblist)
        rvalues = self._parseRvaluefromPDB(pdblist)
        
        # fetchAsymIDs to convert auth_asym_ids into label_asym_id
        auth2label = self.fetchAsymIDs(pdblist)
        
        for pdbid in PDBdata:
            PDBdata[pdbid]['ligand'] = ligands[pdbid]
            PDBdata[pdbid]['initial_release_date'] = rvalues[pdbid]['initial_release_date']
            PDBdata[pdbid]['ls_R_factor_R_free'] = rvalues[pdbid]['ls_R_factor_R_free']
            PDBdata[pdbid]['ls_R_factor_R_work'] = rvalues[pdbid]['ls_R_factor_R_work']
            PDBdata[pdbid]['ls_R_factor_obs'] = rvalues[pdbid]['ls_R_factor_obs']
            
            # fetchAsymIDs to convert auth_asym_ids into label_asym_id
            chains = PDBdata[pdbid]['chains']
            chain_dict = auth2label.get(pdbid, None)
            if chain_dict is None:
                LOGGER.warn(f'fetchAsymIDs: No infor. of {pdbid}')
            else:
                chains = [chain_dict.get(chid, chid) for chid in chains]
            
            # RCSB Data API: polymer_entity_instances 
            # Retrieved info: Sequence Annotations - UNOBSERVED_RESIDUE_XYZ
            pdb_instances = [f'{pdbid}.{chid}' for chid in chains]
            if len(pdb_instances) > 0:
                seq_annot = self._parseSeqAnnotfromPDB(pdb_instances)
            else:
                seq_annot = None
            PDBdata[pdbid]['seq_annot'] = seq_annot
        self._pdbdata = PDBdata
    
    def fetchAsymIDs(self, pdblist):
        """
        Convert auth_asym_ids into label_asym_id. 
        https://www.rcsb.org/docs/general-help/identifiers-in-pdb#:~:text=type%20of%20entity.-,Macromolecular%20Instance%20ID,R%2C%20while%20the%20PDB%20assigned%20ones%20are%20C%20and%20D%20respectively.,-The%20polymer%20sequences
        
        Return: dict
            given auth_asym_ids, dict can be used as a look up table to find label_asym_id.
        """
        query = Query(
            input_type="entries",
            input_ids=pdblist,
            return_data_list=[
                "polymer_entities.polymer_entity_instances.rcsb_id",
                "polymer_entities.rcsb_polymer_entity_container_identifiers.auth_asym_ids",
            ]
        )
        r = query.exec()

        if len(r['data']['entries']) == 0:
            LOGGER.warn(f'fetchAsymIDs: Check input {pdblist}')
            return {}

        auth2label = {}
        for entry in r['data']['entries']:
            pdbid = entry['rcsb_id']
            
            _dict = {}
            for entity in entry['polymer_entities']:
                label_asym_id = entity['polymer_entity_instances'][0]['rcsb_id'].split('.')[-1]
                auth_asym_id  = entity['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'][0]
                _dict[auth_asym_id] = label_asym_id
            
            auth2label[pdbid] = _dict
        return auth2label

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

comma_splitter = re.compile(r'\s*,\s*').split
ns = {'up': 'http://uniprot.org/uniprot'}
    
class UniprotRecord(object):
    """This class provides a wrapper for UniProt data including functions 
    for accessing particular fields and parsing associated PDB entries."""
    def __init__(self, data):
        self._rawdata = data
        self._pdbdata = []
        self._parse()

    def __repr__(self):
        return '<UniprotRecord: %s>'%self.getTitle()

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
        """accession tag"""
        return self.getEntry('accession', index)
    
    def getName(self, index=0):
        """name tag"""
        return self.getEntry('name', index)

    def getProtein(self, index=0):
        """protein tag
        <protein>
            <recommendedName>
                <fullName>Gap junction beta-2 protein</fullName>
            </recommendedName>
            <alternativeName>
                <fullName evidence="57 58">Connexin-26</fullName>
                <shortName evidence="58">Cx26</shortName>
            </alternativeName>
        </protein>
        """
        protein = self.getEntry('protein', index)
        
        try:
            recommend_elem = protein.find('up:recommendedName', ns)
            alternative_elem = protein.find('up:recommendedName', ns)
            
            recommend_name = recommend_elem.find('up:fullName', ns)
            alter_fullname = alternative_elem.find('up:fullName', ns)
            alter_shortname = alternative_elem.find('up:shortName', ns)
            
            recommend_name = recommend_name.text if recommend_name is not None else None
            alter_fullname = alter_fullname.text if alter_fullname is not None else None
            alter_shortname = alter_shortname.text if alter_shortname is not None else None
        except:
            submitted_name = protein.find('up:submittedName/up:fullName', ns)
            submitted_name = submitted_name.text if submitted_name is not None else None
            recommend_name = submitted_name
            alter_fullname = None
            alter_shortname = None
        
        return {
            'recommend_name': recommend_name,
            'alter_fullname': alter_fullname,
            'alter_shortname': alter_shortname            
        }
        
    def getGene(self, index=0):
        """gene tag
        <gene>
            <name type="primary">GJB2</name>
        </gene>
        """
        try:
            gene = self.getEntry('gene', index)
            name_elem = gene.find('up:name[@type="primary"]', ns)
            return name_elem.text if name_elem is not None else None
        
        except Exception as e:
            LOGGER.warn(f'Error while parsing {id}: {e} -> None')
            return None

    def getOrganism(self, index=0):
        """organism tag
        <organism>
            <name type="scientific">Homo sapiens</name>
            <name type="common">Human</name>
            <dbReference type="NCBI Taxonomy" id="9606"/>
            <lineage>
                <taxon>Eukaryota</taxon>
                <taxon>Metazoa</taxon>
                <taxon>Chordata</taxon>
                <taxon>Craniata</taxon>
                <taxon>Vertebrata</taxon>
                <taxon>Euteleostomi</taxon>
                <taxon>Mammalia</taxon>
                <taxon>Eutheria</taxon>
                <taxon>Euarchontoglires</taxon>
                <taxon>Primates</taxon>
                <taxon>Haplorrhini</taxon>
                <taxon>Catarrhini</taxon>
                <taxon>Hominidae</taxon>
                <taxon>Homo</taxon>
            </lineage>
        </organism>
        """
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
        
    def getReference(self, index=0):
        """reference tag"""
        pass
    
    def getComment(self, index=0):
        """conmment tag"""
        pass
        
    def getDBreference(self, index=0):
        """dbReference tag"""
        pass
    
    def getProteinExistence(self, index=0):
        """proteinExistence tag"""
        pass
     
    def getKeyword(self, index=0):
        """keyword tag"""
        pass
    
    def getFeature(self, index=0):
        """feature tag"""
        pass
    
    def getEvidence(self, index=0):
        """evidence tag"""
        pass
    
    def getSequence(self, index=0):
        return self.getEntry('sequence', index)
    
    def getZincFinger(self):
        return self._zinc_finger
    
    def getDNAbinding(self):
        return self._dna_binding
    
    def getActivateSite(self):
        return self._active_site
    
    def getBindingSite(self):
        return self._binding_site
    
    def getSite(self):
        return self._site
    
    def getAlphaFold(self):
        """<dbReference type="AlphaFoldDB" id="Q9UDY8"/>"""
        AlphaFoldDB = None
        for key, value in self._rawdata.items():
            if not key.startswith('dbReference'):
                continue

            if type(value) != list or len(value) != 2:
                continue
            
            # [('type', 'AlphaFoldDB'), ('id', 'Q13509')]
            if value[0][1] == 'AlphaFoldDB':
                AlphaFoldDB = value[1][1]
                break
        return AlphaFoldDB
    
    def getCofactor(self):
        return self._cofactors
    
    def getTitle(self):
        uid = self.getAccession()
        name = self.getName()
        return '%s (%s)'%(uid, name)

    def getEntry(self, item, index=0):
        key = '%s%4d'%(item, index)
        if key in self._rawdata:
            return self._rawdata[key]
        else:
            raise KeyError('%s does not exist in the Uniprot record'%key)

    def _parseDNAbinding(self):
        data = self._rawdata
        dna_binding = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            if value.get('type') != "DNA-binding region":
                continue
            """
            <feature type="DNA-binding region" description="HMG box 1" evidence="4">
                <location>
                <begin position="9"/>
                <end position="79"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            dna_binding.append({
                'description': descp, 
                'begin': begin, 
                'end': end
            })
        self._dna_binding = dna_binding
            
    def _parseZincfinger(self):
        data = self._rawdata
        zinc_finger = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            if value.get('type') != "zinc finger region":
                continue
            """
            <feature type="zinc finger region" description="C2H2-type 1" evidence="1">
                <location>
                <begin position="110"/>
                <end position="133"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            zinc_finger.append({
                'description': descp, 
                'begin': begin, 
                'end': end
            })
        self._zinc_finger = zinc_finger

    def _parseActiveSite(self):
        data = self._rawdata
        active_site = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            
            if value.get('type') != "active site":
                continue
            """
            <feature type="active site" description="Proton donor" evidence="2">
                <location>
                <position position="613"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            pos_elem = value.find('up:location/up:position', ns)
            pos   = int(pos_elem.attrib.get('position')) if pos_elem is not None else None
            active_site.append({
                'description': descp, 
                'position': pos
            })
        self._active_site = active_site
    
    def _parseBindingSite(self):
        data = self._rawdata
        binding_site = []
        for key, value in data.items():
            if not key.startswith('feature'):
                continue
            
            if value.get('type') != "binding site":
                continue
            
            """
            <feature type="binding site" evidence="7 9 22 23 24">
                <location>
                <position position="617"/>
                </location>
                <ligand>
                <name>Zn(2+)</name>
                <dbReference type="ChEBI" id="CHEBI:29105"/>
                </ligand>
            </feature>
            """
            pos_elem = value.find('up:location/up:position', ns)
            pos   = int(pos_elem.attrib.get('position')) if pos_elem is not None else None
            
            ligand_elem = value.find('up:ligand', ns)
            ligand_name = ligand_elem.find('up:name', ns)
            ligand_name = ligand_name.text if ligand_name is not None else None
            ligand_chebi= ligand_elem.find('up:dbReference[@type="ChEBI"]', ns)
            ligand_chebi = ligand_chebi.attrib['id'] if ligand_chebi is not None else None
            binding_site.append({
                'position': pos, 
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
            
            """
            <feature type="site" description="Breakpoint for translocation to form BIRC2-MALT1">
                <location>
                    <begin position="323"/>
                    <end position="324"/>
                </location>
            </feature>
            """
            descp = value.get('description')
            begin_elem = value.find('up:location/up:begin', ns)
            end_elem = value.find('up:location/up:end', ns)
            begin = begin_elem.attrib.get('position') if begin_elem is not None else None
            end = end_elem.attrib.get('position') if end_elem is not None else None
            site.append({
                'description': descp, 
                'begin': begin, 
                'end': end
            })
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
            """
            <comment type="cofactor">
                <cofactor evidence="2">
                    <name>cf_name</name>
                    <dbReference type="ChEBI" id=cf_chebi/>
                </cofactor>
            </comment>
            """
            cf_elem = value.find('up:cofactor', ns)
            # ---
            cf_name = cf_elem.find('up:name', ns)
            cf_chebi= cf_elem.find('up:dbReference[@type="ChEBI"]', ns)
            cf_chebi = cf_chebi.attrib['id'] if cf_chebi is not None else None
            # ---
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
            
            type_list = ['topological domain', 'transmembrane region', 'intramembrane region']
            loc_type  = value.get('type')
            if loc_type not in type_list:
                continue
            """
            <feature type="intramembrane region" evidence="45">
                <location>
                <begin position="2"/>
                <end position="13"/>
                </location>
            </feature>
            """
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

    def _parseSeqAnnotfromPDB(self, pdb_instances):
        """pdb_instances ['1AID.A', '1AID.B']
        """
        query = Query(
            input_type="polymer_entity_instances",
            input_ids=pdb_instances,
            return_data_list=[
                "polymer_entity_instances.rcsb_id",
                "rcsb_polymer_instance_feature.type",
                "rcsb_polymer_instance_feature.feature_positions.beg_seq_id",
                "rcsb_polymer_instance_feature.feature_positions.end_seq_id",
                
                'rcsb_polymer_instance_info.modeled_residue_count',
                'rcsb_polymer_instance_feature_summary.coverage',
                'rcsb_polymer_instance_feature_summary.type',
                'polymer_entity.rcsb_polymer_entity.pdbx_mutation',
                'polymer_entity.entity_poly.rcsb_mutation_count',
            ]
        )
        r = query.exec()

        out = {}
        for entry in r['data']['polymer_entity_instances']:
            instance = entry['rcsb_id']
            instance_feature = entry.get('rcsb_polymer_instance_feature', None)
                
            unobs_res = [ele['feature_positions'] for ele in instance_feature if ele['type'] == 'UNOBSERVED_RESIDUE_XYZ']
            unobs_atom= [ele['feature_positions'] for ele in instance_feature if ele['type'] == 'UNOBSERVED_ATOM_XYZ']
            
            # Turn into a list of residues / atoms
            unobs_res = [
                i for block in unobs_res for r in block
                for i in range(r['beg_seq_id'], r['end_seq_id'] + 1)
            ]
            unobs_atom = [
                i for block in unobs_atom for r in block
                for i in range(r['beg_seq_id'], r['end_seq_id'] + 1)
            ]
            
            # Extract specific types
            coverage_field = entry['rcsb_polymer_instance_feature_summary']    
            unobs_res_cov = next((d['coverage'] for d in coverage_field if d['type'] == 'UNOBSERVED_RESIDUE_XYZ'), 0)
            unobs_atom_cov = next((d['coverage'] for d in coverage_field if d['type'] == 'UNOBSERVED_ATOM_XYZ'), 0)

            coverage = {
                "unobs_res": unobs_res_cov,
                "unobs_atom": unobs_atom_cov
            }
            
            modeled_residue_count = entry['rcsb_polymer_instance_info']['modeled_residue_count']
            pdbx_mutation = entry['polymer_entity']['rcsb_polymer_entity']['pdbx_mutation']
            rcsb_mutation_count = entry['polymer_entity']['entity_poly']['rcsb_mutation_count']
            
            out[instance] = {
                'unobs_res': unobs_res,
                'unobs_atom': unobs_atom,
                'coverage': coverage,
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
            # Return None if no ligand
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
            # Return None if no ligand
            if rcsb_accession_info is None:
                initial_release_date = None
            else:
                initial_release_date = rcsb_accession_info['initial_release_date']
                
            refine = entry.get('refine', None)
            # Return None if no refine
            if refine is None:
                ls_R_factor_R_free = None
                ls_R_factor_R_work = None
                ls_R_factor_obs = None
            else:
                refine_info = refine[0]
                ls_R_factor_R_free = refine_info['ls_R_factor_R_free']
                ls_R_factor_R_work = refine_info['ls_R_factor_R_work']
                ls_R_factor_obs    = refine_info['ls_R_factor_obs']
                
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
            except (KeyError, TypeError) as e:
                continue
            """
            <dbReference type="PDB" id=pdbid>
                <property type="method" value="EM"/>
                <property type="resolution" value=resolution/>
                <property type="chains" value=pdbchains/>
    		</dbReference>
            """
            method = value['method']
            method = value.get('method', None)
            # pdbchains = value['chains'] # e.g. "B/D/F/G/H/I=1-450"
            pdbchains = value.get('chains', []) # e.g. "B/D/F/G/H/I=1-450"
            resolution = value.get('resolution', '1.00 A')
            resolution = float(resolution.split(' ')[0])
            
            # example chain strings: "A=27-139, B=140-150" or "A/B=27-150"
            chains = []
            resrange = None
            try:
                pdbchains = comma_splitter(pdbchains)
                for chain in pdbchains:
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
        if len(pdblist) == 0:
            self._pdbdata = PDBdata
            return
        
        # RCSB Data API: entries 
        # Retrieved info: ligands, released date, Observed Residual factor (R-value obs) 
        ligands = self._parseLigandsfromPDB(pdblist)
        rvalues = self._parseRvaluefromPDB(pdblist)
        
        # fetchAsymIDs to convert auth_asym_ids into label_asym_id
        auth2label = self.fetchAsymIDs(pdblist)
        
        for pdbid in PDBdata:
            PDBdata[pdbid]['ligand'] = ligands[pdbid]
            PDBdata[pdbid]['initial_release_date'] = rvalues[pdbid]['initial_release_date']
            PDBdata[pdbid]['ls_R_factor_R_free'] = rvalues[pdbid]['ls_R_factor_R_free']
            PDBdata[pdbid]['ls_R_factor_R_work'] = rvalues[pdbid]['ls_R_factor_R_work']
            PDBdata[pdbid]['ls_R_factor_obs'] = rvalues[pdbid]['ls_R_factor_obs']
            
            # fetchAsymIDs to convert auth_asym_ids into label_asym_id
            chains = PDBdata[pdbid]['chains']
            chain_dict = auth2label.get(pdbid, None)
            if chain_dict is None:
                LOGGER.warn(f'fetchAsymIDs: No infor. of {pdbid}')
            else:
                chains = [chain_dict.get(chid, chid) for chid in chains]
            
            # RCSB Data API: polymer_entity_instances 
            # Retrieved info: Sequence Annotations - UNOBSERVED_RESIDUE_XYZ
            pdb_instances = [f'{pdbid}.{chid}' for chid in chains]
            if len(pdb_instances) > 0:
                seq_annot = self._parseSeqAnnotfromPDB(pdb_instances)
            else:
                seq_annot = None
            PDBdata[pdbid]['seq_annot'] = seq_annot
        self._pdbdata = PDBdata
    
    def fetchAsymIDs(self, pdblist):
        """
        Convert auth_asym_ids into label_asym_id. 
        https://www.rcsb.org/docs/general-help/identifiers-in-pdb#:~:text=type%20of%20entity.-,Macromolecular%20Instance%20ID,R%2C%20while%20the%20PDB%20assigned%20ones%20are%20C%20and%20D%20respectively.,-The%20polymer%20sequences
        
        Return: dict
            given auth_asym_ids, dict can be used as a look up table to find label_asym_id.
        """
        query = Query(
            input_type="entries",
            input_ids=pdblist,
            return_data_list=[
                "polymer_entities.polymer_entity_instances.rcsb_id",
                "polymer_entities.rcsb_polymer_entity_container_identifiers.auth_asym_ids",
            ]
        )
        r = query.exec()

        if len(r['data']['entries']) == 0:
            LOGGER.warn(f'fetchAsymIDs: Check input {pdblist}')
            return {}

        auth2label = {}
        for entry in r['data']['entries']:
            pdbid = entry['rcsb_id']
            
            _dict = {}
            for entity in entry['polymer_entities']:
                label_asym_id = entity['polymer_entity_instances'][0]['rcsb_id'].split('.')[-1]
                auth_asym_id  = entity['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'][0]
                _dict[auth_asym_id] = label_asym_id
            
            auth2label[pdbid] = _dict
        return auth2label

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
        self._parsePDB()
        LOGGER.report(f'Parsing in %.1fs.', '_parse')
        
