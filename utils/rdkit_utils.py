
import pandas as pd

from base64 import b64encode
import sys
import types
import logging

import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import SDWriter
from rdkit.Chem import rdchem
from rdkit.Chem.Scaffolds import MurckoScaffold
from io import BytesIO
from xml.dom import minidom
from xml.parsers.expat import ExpatError

log = logging.getLogger(__name__)

def LoadSDF(filename, idName='ID', molColName='ROMol', includeFingerprints=False,
            isomericSmiles=True, smilesName=None, embedProps=False, removeHs=True,
            strictParsing=True, chunksize=1000):
    '''Read file in SDF format and return as Pandas data frame.
    If embedProps=True all properties also get embedded in Mol objects in the molecule column.
    If molColName=None molecules would not be present in resulting DataFrame (only properties
    would be read).
    '''
    if isinstance(filename, str):
        if filename.lower()[-3:] == ".gz":
            import gzip
            f = gzip.open(filename, "rb")
        else:
            f = open(filename, 'rb')
            close = f.close
    else:
        f = filename
        close = None  # don't close an open file that was passed in
    
#     RenderImagesInAllDataFrames(images=True)
    
    records = []
    indices = []
    for i, mol in enumerate(
        Chem.ForwardSDMolSupplier(f, sanitize=(molColName is not None), removeHs=removeHs,
                                strictParsing=strictParsing)):
        if mol is None:
            continue
        row = dict((k, mol.GetProp(k)) for k in mol.GetPropNames())
        if molColName is not None and not embedProps:
            for prop in mol.GetPropNames():
                mol.ClearProp(prop)
        if mol.HasProp('_Name'):
            row[idName] = mol.GetProp('_Name')
        if smilesName is not None:
            try:
                row[smilesName] = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)
            except:
                log.warning('No valid smiles could be generated for molecule %s', i)
                row[smilesName] = None
        if molColName is not None and not includeFingerprints:
            row[molColName] = mol
        # elif molColName is not None:
        #     row[molColName] = _MolPlusFingerprint(mol)
        records.append(row)
        indices.append(i)
        
        if len(records) == chunksize:
            yield pd.DataFrame(records, index=indices)
            records = []
            indices = []

    if close is not None:
        close()

    if len(records) > 0:
        assert len(records) < chunksize
        yield pd.DataFrame(records, index=indices)
