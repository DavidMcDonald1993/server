import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

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
from rdkit.Chem.PandasTools import WriteSDF, AddMoleculeColumnToFrame

from io import BytesIO
from xml.dom import minidom
from xml.parsers.expat import ExpatError

from utils.io import read_smiles, smiles_to_sdf

# from standardiser import standardise

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


if __name__ == "__main__":

    smiles_to_sdf("/home/david/Desktop/test_compounds.smi", "/home/david/Desktop/test_compounds.sdf")

    # mol = Chem.MolFromSmiles("Cn1cc(Cn2cnc(-c3cnn(C)c3)c2-c2ccc(C#N)cc2)cn1")
    # print(dict((k, mol.GetProp(k)) for k in mol.GetPropNames()))

    # import glob 

    # smiles_files = glob.iglob(os.path.join("..", "..",
    #     "ppb2", "splits", "*", "test.smi"))

    # sdf_dir = os.path.join("/", "home", "david",
    #     "Desktop", "ppb2_SDF")
    # os.makedirs(sdf_dir, exist_ok=True)

    # for smiles_file in smiles_files:
    #     split_name = smiles_file.split("/")[-2]
    #     sdf_filename = os.path.join(sdf_dir, split_name+".sdf")
    #     smiles_to_sdf(smiles_file, sdf_filename)

        # smiles = read_smiles(smiles_file)

        # smiles["mols"] = smiles["SMILES"].map(lambda smi: standardise_mol(Chem.MolFromSmiles(smi)))

        # smiles["SMILES_standard"] = smiles["mols"].map(Chem.MolToSmiles)

        # print (smiles["SMILES"].isnull().sum())
        # print (smiles["SMILES_standard"].isnull().sum())

        # break

