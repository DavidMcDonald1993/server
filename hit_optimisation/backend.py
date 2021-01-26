import os

import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import numpy as np

import json

from Bio.PDB import PDBParser, PDBIO, PDBList
from Bio.PDB.PDBIO import Select

import scoria

import shutil

from utils.io import process_input_file, BRICS_decompose_smiles_file
from utils.email_utils import send_mail
from utils.users import send_file_to_user, determine_identifier

class ChainSelect(Select):
    
    def __init__(self, chain_id):
        super(ChainSelect, self).__init__()
        self.chain_id = chain_id
    
    def accept_chain(self, chain):
        return chain.id == self.chain_id

def download_pdb_file(
        pdb_id, 
        download_dir):
    print ("downloading PDB file for", pdb_id)

    download_dir = os.path.join(download_dir, 
        "targets", "raw", pdb_id)
    os.makedirs(download_dir, exist_ok=True)
    pdbl = PDBList()
    print ("downloading PDB file to",
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

def remove_residues(pdb_id, 
    pdb_file,
    prepared_dir):
    
    print("preparing PDB file by removing residues")
    prepared_dir = os.path.join(prepared_dir, 
        "targets", "prepared", pdb_id)
    os.makedirs(prepared_dir, exist_ok=True)
    print ("outputting prepared targets to", prepared_dir)

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
        chain_filenames = {}
        for chain in model:
            # print ("writing chain id", chain.id)
            chain_filename = os.path.join(prepared_dir, 
                "{}-chain-{}.pdb".format(pdb_id, chain.id))
            # print ("writing to", chain_filename)
            io.save(chain_filename, ChainSelect(chain.id))

            chain_filenames.update({chain.id: chain_filename})

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

def load_base_settings(
    settings_file=os.path.join("hit_optimisation", "static",
        "hit_optimisation", "base_settings", "base_settings.json")):
    assert os.path.exists(settings_file)
    print ("loading base settings from", settings_file)
    with open(settings_file, "r") as f:
        return json.load(f)

def load_parameter_descriptions(
    settings_file=os.path.join("hit_optimisation", "static",
        "hit_optimisation", "base_settings", "autogrow_parameters.json")):
    assert os.path.exists(settings_file)
    print ("loading base settings from", settings_file)
    with open(settings_file, "r") as f:
        return json.load(f)

def determine_settings(
    chain_filename,
    input_smiles_file,
    bounding_box,
    user_settings,
    output_dir
    ):

    # loading basic settings
    settings = load_base_settings()

    # add receptor to settings
    settings["filename_of_receptor"] = chain_filename

    # add smiles file to setting
    settings["source_compound_file"] = input_smiles_file

    # add bounding box (currently blind docking)
    settings.update(bounding_box)

    # set output dir
    settings["root_output_folder"] = output_dir

    # add user settings
    settings.update(user_settings)

    print ("using settings:")
    for k, v in settings.items():
        assert v is not None, k
        print (k, ":", v)

    settings_filename = os.path.join(output_dir,
        "settings.json")
    save_settings(settings, settings_filename)

    return settings_filename

def save_settings(settings, filename):
    assert filename.endswith(".json")
    print ("writing settings to", filename)

    with open(filename, "w") as f:
        json.dump(settings, f, sort_keys=True, indent=4)


def hit_optimisation(
    user,
    pdb_id, 
    input_file, 
    chain,
    user_settings={},
    compression="zip",
    root_dir="user_files"
    ):
    ''' will run as a process ''' 

    assert len(pdb_id) == 4 

    root_dir = os.path.join(root_dir, 
        "user_id={}".format(user.id), "hit_optimisation")
    os.makedirs(root_dir, exist_ok=True)

    identifier = determine_identifier(input_file)
    identifier += f"-{pdb_id}"

    # process output directory
    output_dir = os.path.join(root_dir, identifier)
    os.makedirs(output_dir, exist_ok=True)
    print ("outputting to directory", output_dir)

    # download file if file is a GET request
    input_smiles_file = process_input_file(input_file, 
        desired_format=".smi", output_dir=output_dir)

    print ("performing hit optimisation for target", pdb_id,
        "with smiles file", input_smiles_file, )

    # perform BRICS decomposition
    # decomposition_file = os.path.join(output_dir,
    #     "BRICS_decomposition.smi")
    input_base, ext = os.path.splitext(input_smiles_file)
    brics_decomposition_out_file = f"{input_base}-BRICS-frags.smi" 
    BRICS_decompose_smiles_file(input_smiles_file, brics_decomposition_out_file)

    # download and process pdb file
    pdb_file = download_pdb_file(pdb_id, 
        download_dir=output_dir)
    assert os.path.exists(pdb_file)

    print ("processing PDB file for AutoDOCK")
    chain_filenames = remove_residues(pdb_id, 
        pdb_file, prepared_dir=output_dir)

    # remove PDB file
    # print ("removing", pdb_file)
    # os.remove(pdb_file)

    print ("using chain", chain)
    chain_filename = chain_filenames[chain]

    # determinine bounding box
    bounding_box = prepare_target_binding_site(chain_filename)

    # determine run settings
    settings_filename = determine_settings(
        chain_filename=chain_filename, 
        # input_smiles_file=input_smiles_file, 
        input_smiles_file=brics_decomposition_out_file,
        bounding_box=bounding_box, 
        user_settings=user_settings, 
        output_dir=output_dir)

    cmd = f'''
    python hit_optimisation/autogrow4/RunAutogrow.py\
        --json "{settings_filename}"
    '''
        # > {output_dir}/autogrow.out\
        # 2> {output_dir}/autogrow.err
    # '''
    
    ret = os.system(cmd)
    # assert ret == 0
    if ret != 0:
        print ("error in run")

    # build zip file containing all targets / run settings / run output
    archive_filename = os.path.join(root_dir, identifier)
    print ("writing archive to", 
        archive_filename + "." + compression)

    shutil.make_archive(archive_filename, 
        compression, output_dir)

    attachment_filename = f"{archive_filename}.{compression}"
    send_file_to_user(user, attachment_filename, 
        subject="NPAIEngine Hit Optimisation Results")

    return 0

class User(object):

    def __init__(self, id, username, email):
        self.id = id
        self.username = username 
        self.email = email

if __name__ == "__main__":

    user_id = 1
    username = "David"
    email = "davemcdonald93@gmail.com"
    pdb_id = "4DQY"
    smiles_file = "./all_ligands.smi"
    chain = "C"

    user = User(user_id, username, email)

    hit_optimisation(
        user, 
        pdb_id, 
        smiles_file, 
        chain,
        user_settings={"num_generations": 1}
    )