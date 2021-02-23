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

from rdkit.Chem.PandasTools import LoadSDF, WriteSDF
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.PandasTools import WriteSDF, AddMoleculeColumnToFrame
from rdkit.Chem.BRICS import BRICSDecompose

import re

def sanitise_filename(filename):
    return filename.replace(" ", "")

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

    try:
        AllChem.EmbedMolecule(my_mol_with_H,useRandomCoords=False)
        AllChem.MMFFOptimizeMolecule(my_mol_with_H)
    except ValueError:
        AllChem.EmbedMolecule(my_mol_with_H,useRandomCoords=True)
        AllChem.MMFFOptimizeMolecule(my_mol_with_H)


    embedded_mol = Chem.RemoveHs(my_mol_with_H)
    return embedded_mol

def smiles_to_sdf(
    smiles_filename, 
    sdf_filename,
    standardise=True,
    embed=False):
    print ("converting smiles from", smiles_filename, 
        "to SDF file", sdf_filename)
    smiles_df = read_smiles(smiles_filename, )

    print ("num smiles:", smiles_df.shape[0])

    AddMoleculeColumnToFrame(smiles_df, 'SMILES', 'Molecule')
    molColName = "Molecule"

    if standardise:
        print ("standardising SMILES")
        smiles_df["MoleculeStandard"] = smiles_df["SMILES"].map(standardise_smi, na_action="ignore")
        smiles_df["SMILESStandard"] = smiles_df["MoleculeStandard"].map(Chem.MolToSmiles, na_action="ignore")
        molColName = "MoleculeStandard"

    if embed:
        print ("embedding SMILES into 3D")
        smiles_df["MoleculeEmbedded"] = smiles_df["SMILES"].map(embed_2D_mol_in_3D, na_action="ignore")
        molColName = "MoleculeEmbedded"

    smiles_df = smiles_df.loc[~pd.isnull(smiles_df[molColName])] # drop missing values
    print ("num SMILES remaining:", smiles_df.shape[0])
    WriteSDF(smiles_df, sdf_filename, molColName=molColName,
        idName="RowID", properties=list(smiles_df.columns))

def copy_file(input_file, output_dir):
    print ("input file is a local file")
    input_file_sanitised = sanitise_filename(input_file)
    input_file_type = os.path.splitext(input_file_sanitised)[1]
    assert input_file_type in valid_input_file_types
    temp_file = os.path.join(output_dir,
        os.path.basename(input_file_sanitised)) # no uploaded file
    print ("copying file to", temp_file)
    shutil.copyfile(input_file, temp_file)
    return temp_file

def download_from_client(input_file, output_dir):
    assert hasattr(input_file, "name")
    input_filename = input_file.name
    print ("NAME", input_filename)
    input_filename = sanitise_filename(input_filename)
    print ("input file has been received from client -- downloading")

    # write compounds to local directory of server
    temp_file = os.path.join(output_dir,
        input_filename)
    print ("downloading file to", temp_file)
    with open(temp_file, "wb+") as out_file:
        for chunk in input_file.chunks():
            out_file.write(chunk)
    return temp_file


def replace_missing_names(names):
    missing = 1
    changed = False
    for i in range(len(names)):
        if names[i] == "":
            changed = True
            names[i] = f"unknown_compound_{i}"
    return changed, names

def get_compound_names(filename):
    print ("getting compound names from", filename)
    if filename.endswith(".sdf"):
        sdf_id_col = "ID"
        sdf = LoadSDF(filename, idName=sdf_id_col)
        names = list(sdf[sdf_id_col].values)
        changed, names = replace_missing_names(names)
        if changed:
            sdf[sdf_id_col] = names
            print ("names have changed, writing SDF to file", filename)
            WriteSDF(sdf, out=filename, idName=sdf_id_col)
    else: # assume smiles file
        smiles = read_smiles(filename, 
            filter_valid=True, return_series=True)
        names = list(smiles.index)
        _, names = replace_missing_names(names)
        # if changed:
            
        print ("names have changed, writing smiles to file", filename)
        write_smiles([(name, smi) 
            for name, smi in zip(names, smiles.values)], 
        filename)

    return names

def convert_file(
    input_file, 
    desired_format,
    output_dir, 
    valid_input_file_types=(".smi", ".txt", ".sml", ".sdf"),
    ):

    input_file_name, input_file_type = os.path.splitext(input_file)
    assert input_file_type in valid_input_file_types

    # convert if necessary
    if input_file_type != desired_format:
        print ("CONVERSION NECESSARY")
        if desired_format == ".smi":
            print ("COVERTING TO smi")
            if input_file_type == ".sdf":
                # convert SDF to smiles
                print ("converting SDF to SMILES")
                print ("SDF filename:", temp_file)
                sdf_df = LoadSDF(input_file, smilesName="SMILES")
                # write smiles
                smiles = [(row["ID"], row["SMILES"])
                    for _, row in sdf_df.iterrows()]
                # write smiles to temp_file
                input_file = input_file_name + desired_format
                write_smiles(smiles, temp_file)
            elif input_file_type in {".txt", ".sml"}:
                # rename .txt smiles format to .smi
                print ("INPUT FILE TYPE:", input_file_type)
                os.rename(input_file, input_file_name + desired_format)
                input_file = input_file_name + desired_format
            else:
                raise NotImplementedError

        elif desired_format == ".sdf":
            if input_file_type in {".smi", ".sml", ".txt"}:
                # convert from SMILES to SDF
                print ("converting SMILES to SDF")
                smiles_filename = input_file
                print ("SMILES filename:", smiles_filename)
                input_file = input_file_name + desired_format
                smiles_to_sdf(smiles_filename, input_file,
                    standardise=True, embed=False)
            else:
                raise NotImplementedError
        
        else:
            raise NotImplementedError #conversion not yet implemented
    else:
        print ("no conversion necessary")
    
    assert input_file.endswith(desired_format)

    return input_file

def BRICS_decompose_smiles(smi, minFragmentSize=3):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return mol
    return BRICSDecompose(mol, minFragmentSize=minFragmentSize)

def BRICS_decompose_smiles_file(smiles_file, out_file, keep_original=True):
    print ("performing BRICS decomposition on smiles file", smiles_file,
        "outputting to", out_file)
    smiles = read_smiles(smiles_file, return_series=True, )

    expr = re.compile(r'\[[0-9]+\*\]')
    empty_brackets = re.compile(r"\(\)")

    decomposed_smiles = dict()

    for compound, smi in smiles.items():
        if keep_original and smi not in decomposed_smiles:
            decomposed_smiles[smi] = compound
        for i, fragment in enumerate(map( lambda s: 
            empty_brackets.sub("", expr.sub("", s)), BRICS_decompose_smiles(smi))):
            # decomposed_smiles.append((f"{compound}_{i+1}", fragment))
            if fragment not in decomposed_smiles:
                decomposed_smiles[fragment] = f"{compound}_BRICS_frag_{i+1}"

    decomposed_smiles = [(compound, smi) for smi, compound in decomposed_smiles.items()]

    write_smiles(decomposed_smiles, out_file)

if __name__ == "__main__":
    
    input_file = "/home/david/Desktop/Book1.smi"
    desired_format = ".sdf"
    output_dir = "."
    # out_file = "decomposed_test.smi"

    compound_names = get_compound_names(input_file)
    print (compound_names)

    # processed_file = process_input_file(input_file, desired_format, output_dir)

    # BRICS_decompose_smiles_file(input_file, out_file)

    # print (processed_file)

    # smi = "Cn1cc(Cn2cnc(-c3cnn(C)c3)c2-c2ccc(C#N)cc2)cn1"
    # embedded_mol = embed_2D_mol_in_3D(smi)

    # print (Chem.MolToSmiles(embedded_mol))


    # parameter_info = load_json("autogrow_parameters.json")

    # print (parameter_info)