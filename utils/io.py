import json 

def load_json(json_filename):
    print ("loading json from", json_filename)
    with open(json_filename,"r") as f:
        return json.load(f)

def write_smiles(smiles, smiles_filename):
    assert isinstance(smiles, list) # list of (compound_id, smiles) tuples
    print ("writing", len(smiles), "smiles to", smiles_filename)
    with open(smiles_filename, "w") as f:
        for compound_id, smile in smiles:
            f.write(f"{compound_id}\t{smile}\n")

