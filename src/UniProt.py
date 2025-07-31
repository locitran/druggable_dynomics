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

from .utils.logger import LOGGER

comma_splitter = re.compile(r'\s*,\s*').split

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

    def getSequence(self, index=0):
        return self.getEntry('sequence', index)
    
    def getAccession(self, index=0):
        return self.getEntry('accession', index)
    
    def getName(self, index=0):
        return self.getEntry('name', index)

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
            pdbchains = value['chains']
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
                'method': value['method'],
                'resolution': resolution,
                'chains': chains,
            }

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
