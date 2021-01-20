import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import os

import json 

import pandas as pd

from scipy import sparse as sp

import shutil

from standardiser import standardise

from rdkit.Chem.PandasTools import LoadSDF
from rdkit import Chem
from rdkit.Chem.PandasTools import WriteSDF, AddMoleculeColumnToFrame

def valid_smiles(smi):
    assert smi is not None
    try:
        return Chem.MolFromSmiles(smi) is not None
    except TypeError:
        assert False, smi
        return False

def load_json(json_filename):
    print ("loading json from", json_filename)
    with open(json_filename,"r") as f:
        return json.load(f)

def write_json(data, json_filename):
    print ("writing json to", json_filename)
    with open(json_filename, "w") as f:
        json.dump(data, f, indent=4)

def write_smiles(smiles, smiles_filename):
    assert isinstance(smiles, list) # list of (compound_id, smiles) tuples
    print ("writing", len(smiles), "smiles to", smiles_filename)
    with open(smiles_filename, "w") as f:
        for compound_id, smile in smiles:
            f.write(f"{smile}\t{compound_id}\n")

def read_smiles(smiles_filename, filter_valid=True, return_series=False, smiles_col="SMILES"):
    print ("reading smiles from", smiles_filename)
    assert os.path.exists(smiles_filename)
    smiles_df = pd.read_csv(smiles_filename, 
        names=[smiles_col, "compound"],
        sep="\t", header=None)
    smiles_df = smiles_df.loc[~pd.isnull(smiles_df[smiles_col])]
    if filter_valid:
        print ("removing invalid smiles")
        smiles_df = smiles_df.loc[smiles_df[smiles_col].map(valid_smiles)]
    smiles_df = smiles_df.set_index("compound", drop=True)
    if return_series:
        smiles_df = smiles_df[smiles_col]
    return smiles_df

def load_labels(labels_filename):
    print ("loading labels from", labels_filename)
    Y = sp.load_npz(labels_filename)
    print ("labels shape is", Y.shape)
    return Y # sparse format

def standardise_smi(smi, return_smiles=False):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        if return_smiles:
            return smi 
        else:
            return mol
    try:
        mol = standardise.run(mol)
    except standardise.StandardiseException as e:
        print (e)
        pass
    if return_smiles:
        return Chem.MolToSmiles(mol)
    else:
        return mol

def embed_2D_mol_in_3D(smi):
    assert smi is not None
    mol = Chem.MolFromSmiles(smi)
    my_mol_with_H = Chem.AddHs(mol)

    AllChem.EmbedMolecule(my_mol_with_H)
    AllChem.MMFFOptimizeMolecule(my_mol_with_H)

    embedded_mol = Chem.RemoveHs(my_mol_with_H)
    return embedded_mol

def smiles_to_sdf(
    smiles_filename, 
    sdf_filename):
    print ("converting smiles from", smiles_filename, 
        "to SDF file", sdf_filename)
    smiles_df = read_smiles(smiles_filename, )
    smiles_df["Molecule"] = smiles_df["SMILES"].map(embed_2D_mol_in_3D)
    smiles_df["MoleculeStandard"] = smiles_df["SMILES"].map(standardise_smi)
    smiles_df["clean_SMILES"] = smiles_df["MoleculeStandard"].map(Chem.MolToSmiles, na_action="ignore")
    smiles_df = smiles_df.loc[~pd.isnull(smiles_df["Molecule"])] # drop missing values
    # AddMoleculeColumnToFrame(smiles_df, 'SMILES', 'Molecule')
    WriteSDF(smiles_df, sdf_filename, molColName="Molecule",
        idName="RowID", properties=list(smiles_df.columns))

def process_input_file(
    input_file, 
    desired_format,
    output_dir, 
    valid_input_file_types=(".smi", ".txt", ".sml", ".sdf")):
    '''
    process input file from client
    '''
    assert desired_format in valid_input_file_types 
    if not isinstance(input_file, str):
        assert hasattr(input_file, "name")
        print ("input file has been recieved from client -- downloading")
        input_file_type = os.path.splitext(input_file.name)[1]
        assert input_file_type in valid_input_file_types
        # write compounds to server local directory
        temp_file = os.path.join(output_dir,
            input_file.name)
        with open(temp_file, "wb+") as out_file:
            for chunk in input_file.chunks():
                out_file.write(chunk)
    else:
        print ("input file is a local file")
        input_file_type = os.path.splitext(input_file)[1]
        assert input_file_type in valid_input_file_types
        temp_file = os.path.join(output_dir,
            os.path.basename(input_file)) # no uploaded file
        shutil.copyfile(input_file, temp_file)
    
    temp_file_name, _ = os.path.splitext(temp_file)

    # convert if necessary
    if input_file_type != desired_format:
        if desired_format == ".smi":
            if input_file_type == ".sdf":
                # convert SDF to smiles
                print ("converting SDF to SMILES")
                print ("SDF filename:", temp_file)
                sdf_df = LoadSDF(temp_file, smilesName="SMILES")
                # write smiles
                smiles = [(row["ID"], row["SMILES"])
                    for _, row in sdf_df.iterrows()]
                # write smiles to temp_file
                temp_file = temp_file_name + desired_format
                write_smiles(smiles, temp_file)
            elif input_file_type in {".txt", ".sml"}:
                # rename .txt smiles format to .smi
                os.rename(temp_file, temp_file_name + desired_format)
            else:
                raise NotImplementedError

        elif desired_format == ".sdf":
            if input_file_type in {".smi", ".sml", ".txt"}:
                # convert from SMILES to SDF
                print ("converting SMILES to SDF")
                smiles_filename = temp_file
                print ("SMILES filename:", smiles_filename)
                temp_file = temp_file_name + desired_format
                smiles_to_sdf(smiles_filename, temp_file)
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError #conversion not yet implemented
    else:
        print ("no conversion necesary")

    
    assert temp_file.endswith(desired_format)
    
    return temp_file

if __name__ == "__main__":
    
    input_file = "temp.smi"
    desired_format = ".sdf"
    output_dir = "/home/david/Desktop"
    
    processed_file = process_input_file(input_file, desired_format, output_dir)

    print (processed_file)