
import os
import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from django.utils.encoding import smart_str

from datetime import datetime

from pymongo import MongoClient

import numpy as np
import pandas as pd

import shutil

from rdkit.Chem.PandasTools import LoadSDF

# from utils.email_utils import send_mail
# from utils.mongodb_utils import connect_to_mongodb
from utils.pass_utils import (remove_invalid_characters, parse_pass_spectra, 
    get_all_targets, determine_targets)
from utils.enrichment_utils import perform_enrichment_analysis
from utils.io import process_input_file, write_json
from utils.genenames_utils import targets_to_uniprot_ids
from utils.mysql_utils import get_uniprots_for_targets
# from utils.rdkit_utils import LoadSDF
from utils.users import send_file_to_user

# def write_PASS_hits_to_db(
#     pass_file,
#     collection="PASS_hits"):
#     assert isinstance(pass_file, str)
#     assert pass_file.endswith(".sdf")
#     print ("reading PASS hits from", pass_file,
#         "and extracting activity data")

#     pass_activities = LoadSDF(pass_file, 
#         smilesName="SMILES", molColName=None)

#     db = connect_to_mongodb()

#     print ("writing PASS hits to collection:", collection)

#     hit_collection = db[collection]

#     records = []

#     print ("determining hits")

#     for _, row in pass_activities.iterrows():
        
#         entry = {}
#         for col in row.index:
#             if pd.isnull(row[col]):
#                 continue
#             if col == "PASS_ACTIVITY_SPECTRUM":
#                 value = list(
#                     map(parse_pass_spectra, 
#                         map(remove_invalid_characters,
#                             row[col].split("\n"))))
#             else:
#                 value = row[col]
#             entry.update({col: value})
#         entry.update({"time": str(datetime.now())})
#         records.append(entry)

#     print ("inserting", len(records), "records")
#     hit_collection.insert_many(records)

#     db.client.close()

def perform_enrichment_on_PASS_file(
    pass_out_file, 
    output_dir,
    threshold=500,
    ):

    output_dir = os.path.join(output_dir, "enrichment")
    os.makedirs(output_dir, exist_ok=True)

    print ("perfoming enrichment analysis on predicted PASS file",
        "to directory", output_dir)

    # determine active targets from PASS-predicted SDF file
    active_targets = determine_targets(pass_out_file, threshold=threshold)
    active_targets_filename = os.path.join(output_dir,
        f"active_targets_threshold={threshold}.json")
    write_json(active_targets, active_targets_filename)

    unique_target_names = {target 
        for compound, targets in active_targets.items()
        for target in targets}
    
    # unique_uniprots = targets_to_uniprot_ids(unique_target_names)
    targets_to_uniprot = get_uniprots_for_targets(unique_target_names) # SQL query
    targets_to_uniprot_filename = os.path.join(output_dir, "targets_to_uniprot.csv")
    print ("writing targets to uniprot to", targets_to_uniprot_filename)
    targets_to_uniprot = pd.DataFrame(targets_to_uniprot, 
        columns=["target", "uniprot_ACC", "association_score"])
    targets_to_uniprot.to_csv(targets_to_uniprot_filename)

    unique_uniprots = set(targets_to_uniprot["uniprot_ACC"])
    unique_uniprots_filename = os.path.join(output_dir, 
        "unique_uniprot_ACCs.txt")
    print ("writing unique uniprots to", unique_uniprots_filename)
    with open(unique_uniprots_filename, "w") as f:
        f.write("\n".join(unique_uniprots))

    # declare filenames to output enrichment
    output_csv_filename = os.path.join(output_dir, 
        "enrichment.csv")
    found_filename = os.path.join(output_dir,
        "found.txt")
    not_found_filename = os.path.join(output_dir,
        "not_found.txt")
    pdf_filename = os.path.join(output_dir,
        "enrichment_summary.pdf")

    perform_enrichment_analysis(
        unique_uniprots_filename,
        output_csv_filename,
        found_filename,
        not_found_filename,
        pdf_filename)

    return 0

def determine_identifier(receiver_address, input_file):
    if not isinstance(input_file, str):
        assert hasattr(input_file, "name")
        input_file = input_file.name 
    # assert smiles_file.endswith(".smi")
    return "{}-{}".format(receiver_address, 
        os.path.splitext(os.path.basename(input_file))[0])

def pass_predict(
    user,
    # user_email,
    input_file, 
    compression="zip",
    static_dir="pass_app/static/pass_app",
    output_dir="files",
    archive_dir="archives",
    threshold=500,
    enrichment=True):

    identifier = determine_identifier(user.email, input_file)

    output_dir = os.path.join(static_dir, output_dir, identifier)
    os.makedirs(output_dir, exist_ok=True)

    archive_dir = os.path.join(static_dir, archive_dir)
    os.makedirs(archive_dir, exist_ok=True)

    input_file = process_input_file(input_file, 
        desired_format=".sdf", output_dir=output_dir)

    base_name, extension = os.path.splitext(input_file)
    assert extension == ".sdf"

    pass_out_file = base_name +"-PASS-out.sdf"

    cmd = "PASS2019toSDF.exe {} {}".format(input_file, pass_out_file)
    print ("executing command:", cmd)

    ret = os.system(cmd)
    assert ret == 0
    assert os.path.exists(pass_out_file)

    # add enrichment to output directory
    if enrichment:
        ret = perform_enrichment_on_PASS_file(pass_out_file,
            output_dir=output_dir, threshold=threshold)

    # build zip file containing all targets / run settings / run output
    archive_filename = os.path.join(archive_dir,
        identifier)
    print ("writing archive to", 
        archive_filename + "." + compression)

    shutil.make_archive(archive_filename, 
        compression, output_dir)

    attachment_filename = f"{archive_filename}.{compression}"
    send_file_to_user(user,attachment_filename)

    return 0

if __name__ == "__main__":
    name = "david"
    email = "davemcdonald93@gmail.com"
    input_file = "/home/david/Desktop/test-PASS-out.sdf"

    # determine_targets(input_file)
    perform_enrichment_on_PASS_file(input_file, output_dir="/home/david/Desktop", threshold=0)



