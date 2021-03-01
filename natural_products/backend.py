import os
import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import re

import numpy as np
import pandas as pd

from datetime import datetime

from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

import urllib.parse as urlparse

from utils.io import load_json, write_json, write_smiles
from utils.mongodb_utils import connect_to_mongodb
from utils.mysql_utils import mysql_query, mysql_create_table, mysql_insert_many, connect_to_mysqldb, sanitise_names
from utils.pass_utils import get_categories, get_all_targets, get_targets_for_category, get_all_compounds
from utils.enrichment_utils import perform_enrichment_analysis


def get_coconut_compound_info_from_mongo( # use mongo
    compound_id, 
    projection={"_id": 0},
    compound_info_collection="uniqueNaturalProduct",
    filter_pa_pi=True,
    ):

    print ("querying COCONUT database for info",
        "about compound with compound id", compound_id)

    db = connect_to_mongodb()

    coconut_collection = db[compound_info_collection]

    compound_info = coconut_collection.find_one(
        {"coconut_id": compound_id},
        projection=projection) # get all info from mongo

    return compound_info


def draw_molecule(
    compound_id,
    smiles, 
    static_dir="static",
    output_dir="compound_images",
    ):
    output_dir = os.path.join(output_dir, f"{compound_id//1024}")
    os.makedirs(os.path.join(static_dir, output_dir), exist_ok=True)

    img_filename = os.path.join(output_dir, f"{compound_id}.png")
    img_full_path = os.path.join(static_dir, img_filename)
    if os.path.exists(img_full_path):
        print (img_full_path, "already exists")
        return img_filename
   
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        img = MolToImage(mol)
        img.save(img_full_path)
        print ("drawimg molecule to", img_full_path)
        return img_filename
    else:
        return None

def write_records_to_file(
    user_id,
    targets, 
    threshold,
    records,
    root_dir="user_files",
    ):
    assert isinstance(records, pd.DataFrame)
    if "image" in records.columns:
        del records["image"]

    output_dir = os.path.join(root_dir, f"user_id={user_id}", "hits")
    os.makedirs(output_dir, exist_ok=True)

    targets = ",".join(map(lambda s: re.sub(r"( |/)", "_", s), targets))
    # thresholds= ",".join(map(str, thresholds))

    records_filename = os.path.join(output_dir,
        f'''targets={targets}-threshold={threshold}-results.csv''')
    print ("writing records to", records_filename)
    records.to_csv(records_filename)

    return records_filename

def write_smiles_to_file(
    user_id,
    targets,
    threshold,
    smiles,
    root_dir="user_files",
    ):  

    output_dir = os.path.join(root_dir,f"user_id={user_id}", "smiles")
    os.makedirs(output_dir, exist_ok=True)

    targets = ",".join(map(lambda s: re.sub(r"( |/)", "_", s), targets))
    # thresholds= ",".join(map(str, thresholds))

    smiles_filename = os.path.join(output_dir,
        f'''targets={targets}-threshold={threshold}-hits.smi''')
    print ("writing smiles to", smiles_filename)

    write_smiles(smiles, smiles_filename)

    return smiles_filename
    
if __name__ == "__main__":

    # from timeit import default_timer

    query = '''
    SELECT coconut_id, clean_smiles
    FROM compounds
    '''
    records = mysql_query(query)

    # for _id, smiles in records:
    #     draw_molecule(_id, smiles)
    write_smiles(records, "coconut_smiles_clean.smi")