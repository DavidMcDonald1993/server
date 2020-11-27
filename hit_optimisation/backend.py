
import os

import numpy as np

import json

from Bio.PDB import PDBParser, PDBIO, PDBList
from Bio.PDB.PDBIO import Select

import scoria

class ChainSelect(Select):
    
    def __init__(self, chain_id):
        super(ChainSelect, self).__init__()
        self.chain_id = chain_id
    
    def accept_chain(self, chain):
        return chain.id == self.chain_id

def download_pdb_file(pdb_id, 
    download_dir=os.path.join("hit_optimisation", 
        "static", "hit_optimisation",
        "targets", "raw")):
    print ("downloading PDB file for", pdb_id)

    download_dir = os.path.join(download_dir, pdb_id)
    os.makedirs(download_dir, exist_ok=True)
    pdbl = PDBList()
    print ("outputting PDB file to",
        download_dir)
    return pdbl.retrieve_pdb_file(pdb_id,
        file_format="pdb", 
        pdir=download_dir)

# def valid_pdb(pdb_file):
#     if not os.path.exists(pdb_file):
#         return False
#     try:
#         parser = PDBParser()
#         parser.get_structure("test", 
#             pdb_file)
#         return True 
#     except Exception as e:
#         return False

def remove_residues(pdb_id, pdb_file,
    prepared_dir=os.path.join("hit_optimisation", 
        "static", "hit_optimisation",
        "targets", "prepared")):
    
    print("preparing PDB file by removing residues")
    prepared_dir = os.path.join(prepared_dir, pdb_id)
    os.makedirs(prepared_dir, exist_ok=True)
    print ("outputting to", prepared_dir)

    parser = PDBParser()
    structure = parser.get_structure(pdb_id, 
        pdb_file)

    if len(structure) == 1: # TODO

        model = structure[0]
        residue_to_remove = []
        chain_to_remove = []

        for chain in model:
            # print ("processing chain id", chain.id)
            for residue in chain:
                if residue.id[0] != ' ':
                    residue_to_remove.append((chain.id, residue.id))
            if len(chain) == 0:
                chain_to_remove.append(chain.id)

        # write for each chain
        io = PDBIO()
        io.set_structure(structure)
        chain_filenames = []
        for chain in model:
            # print ("writing chain id", chain.id)
            chain_filename = os.path.join(prepared_dir, 
                "{}.pdb".format(chain.id))
            # print ("writing to", chain_filename)
            io.save(chain_filename, ChainSelect(chain.id))

            chain_filenames.append((chain.id, chain_filename))

    return chain_filenames
    
def prepare_target_binding_site(pdb_file):
    print ("determining binding site for pdb file",
        pdb_file)
    print ("preparing bounding box for blind docking")

    # create a scoria mol object from the protein pdb file
    mol = scoria.Molecule(pdb_file)

    center_x, center_y, center_z = identify_centre(mol)

    size_x, size_y, size_z = get_bounding_box_size(mol)

    bounding_box = {
        "center_x": center_x,
        "center_y": center_y,
        "center_z": center_z,
        "size_x": size_x,
        "size_y": size_y,
        "size_z": size_z,
    }

    return bounding_box

def identify_centre(mol):
    print ("identifying mol centre of mass")
    return mol.get_center_of_mass().data

def get_bounding_box_size(mol, allowance=2.):
    print ("determining bounding box of mol")
    bounding_box = mol.get_bounding_box()
    return np.ceil(np.abs(bounding_box[0] - bounding_box[1])) + allowance 

def load_settings(
    settings_file=os.path.join("hit_optimisation", "static",
        "hit_optimisation", "settings.json")):
    assert os.path.exists(settings_file)
    print ("loading base settings from", settings_file)

    with open(settings_file, "r") as f:
        return json.load(f)

def hit_optimisation(
    target, 
    smiles_file, 
    num_generations=5,
    output_dir=os.path.join("hit_optimisation",
        "static", "hit_optimisation", "output")):

    assert len(target) == 4 # PDB_ID

    print ("performing hit optimisation for target", target,
        "with smiles file", smiles_file, 
        "using", num_generations, "generations")

    if not isinstance(smiles_file, str):

        assert smiles_file.name.endswith(".smi")

        # print ("performing hit optimisation for target", target,
        #     "using smiles file", smiles_file.name)

        # write compounds to server local directory
        temp_smiles_file = "temp.smi"
        with open(temp_smiles_file, "wb+") as out_file:
            for chunk in smiles_file.chunks():
                out_file.write(chunk)
    else:
        temp_smiles_file = smiles_file # no request

    # process output directory
    output_dir = os.path.join(output_dir, 
        "{}-{}".format(target, temp_smiles_file))
    os.makedirs(output_dir, exist_ok=True)
    print ("outputting to directory", output_dir)
    
    # loading basic settings
    settings = load_settings()

    # download and process pdb file
    pdb_file = download_pdb_file(target)
    assert os.path.exists(pdb_file)

    print ("processing PDB file for AutoDOCK")
    chain_filenames = remove_residues(pdb_id, pdb_file)

    # remove PDB file
    print ("removing", pdb_file)
    os.remove(pdb_file)

    chosen_chain, chain_filename = chain_filenames[0]
    print ("using chain", chosen_chain)

    bounding_box = prepare_target_binding_site(chain_filename)

    # add receptor to settings
    settings["filename_of_receptor"] = chain_filename

    # add smiles file to setting
    settings["source_compound_file"] = temp_smiles_file

    # add bounding box (currently blind docking)
    settings.update(bounding_box)

    # add number of generations
    settings["num_generations"] = num_generations

    # set output dir
    settings["root_output_folder"] = output_dir

    for k, v in settings.items():
        assert v is not None, k

    # TODO finish

    #

if __name__ == "__main__":
    pdb_id = "5FRE"
    smiles_file = "temp.smi"

    # pdb_file = download_pdb_file(pdb_id)
    # assert os.path.exists(pdb_file)

    # chain_filenames = remove_residues(pdb_id, pdb_file)

    # for chain_id, chain_filename in chain_filenames:

    #     bounding_box = prepare_target_binding_site(chain_filename)

    #     print (chain_id, chain_filename)
    #     print (bounding_box)
    #     print ()

    # settings = load_settings()

    # print (settings)

    hit_optimisation(pdb_id, smiles_file)