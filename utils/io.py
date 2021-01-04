import os

import json 

import pandas as pd

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

def read_smiles(smiles_filename):
    print ("reading smiles from", smiles_filename)
    assert os.path.exists(smiles_filename)
    smiles_df = pd.read_csv(smiles_filename, 
        names=["SMILES", "compound"],
        sep="\t", header=None)
    return smiles_df.set_index("compound", drop=True)