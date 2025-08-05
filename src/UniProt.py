# -*- coding: utf-8 -*-
"""This module defines a class and relative functions for mapping Uniprot
sequences to PDB and Pfam databases."""

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
from rcsbapi.data import DataQuery

from .utils.logger import LOGGER

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
        """protein tag"""
        protein = self.getEntry('protein', index)
        protein = protein.find('up:recommendedName/up:fullName', ns)
        if protein is not None:
            return protein.text

    def getGene(self, index=0):
        """gene tag"""
        pass

    def getOrganism(self, index=0):
        """organism tag"""
        pass

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
    
    def _parseCofactor(self):
        data = self._rawdata
        cofactors = []
        for key, value in data.items():
            if not key.startswith('comment'):
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
    
    def _parseLigand(self, pdblist):
        """
        Fetch data from RCSB graphQL using Data API
        https://rcsbapi.readthedocs.io/en/latest/data_api/quickstart.html
        
        # Available return_data_list fields: 
        # https://data.rcsb.org/data-attributes.html 
        """
        
        query = DataQuery(
            input_type="entries",
            input_ids=pdblist,
            return_data_list=[
                "nonpolymer_entities.pdbx_entity_nonpoly.comp_id",
            ]
        )
        r = query.exec()
        
        ligands = {}
        for entry in r['data']['entries']:
            pdbid = entry['rcsb_id']
            ligand_list = entry.get('nonpolymer_entities', None)
            # Return None if no ligand
            if ligand_list is None:
                ligands[pdbid] = None
                continue

            ligand = [
                entity['pdbx_entity_nonpoly']['comp_id'] \
                for entity in ligand_list
            ]
            ligands[pdbid] = ligand
        return ligands
    
    def _parse(self):
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
            pdbchains = value['chains'] # e.g. "B/D/F/G/H/I=1-450"
            resolution = value.get('resolution', '1.00 A')
            resolution = float(resolution.split(' ')[0])
            
            # example chain strings: "A=27-139, B=140-150" or "A/B=27-150"
            chains = []
            pdbchains = comma_splitter(pdbchains)
            for chain in pdbchains:
                chids, resrange = chain.split('=')
                chids = [chid.strip() for chid in chids.split('/')]
                for chid in chids:
                    chains.append(chid)
                    
            PDBdata[pdbid] = {
                'method': method,
                'resolution': resolution,
                'chains': chains,
                'resrange': resrange,
            }

        pdblist = list(PDBdata.keys())
        if len(pdblist) != 0:
            ligands = self._parseLigand(pdblist)
            for pdbid in PDBdata:
                PDBdata[pdbid]['ligand'] = ligands[pdbid]

        self._parseActiveSite()
        self._parseBindingSite()
        self._parseCofactor()
        self._parseDNAbinding()
        self._parseZincfinger()
        
        self._pdbdata = PDBdata
        
    
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
