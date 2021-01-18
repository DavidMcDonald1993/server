
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
from utils.io import process_input_file, write_json, read_smiles, load_json
from utils.genenames_utils import targets_to_uniprot_ids
from utils.queries import get_uniprots_for_targets
# from utils.rdkit_utils import LoadSDF
from utils.users import send_file_to_user, determine_identifier
from utils.ppb2_utils import load_model

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

def perform_predicton_with_novel_classifier(
    # smiles_file,
    smiles,
    model_filename="models/morg3-xgc.pkl.gz",
    n_proc=6,
    max_confidence=1000):
    '''
    Predict with novel classifier
    '''
    # assert smiles_file.endswith(".smi")
    # # read (and filter smiles)
    # smiles = read_smiles(
    #     smiles_file,
    #     filter_valid=True, 
    #     return_series=True)

    model = load_model(model_filename)
    if hasattr(model, "n_proc"):
        model.set_n_proc(n_proc)
    
    # make prediction using pretrained model
    predictions = model.predict_proba(smiles)

    # rescale by max confidence
    predictions = (predictions * max_confidence / predictions.max(axis=0, keepdims=True)).astype(int)

    id_to_db_id = load_json("id_to_db_id.json")

    return pd.DataFrame(predictions, 
        index=smiles.index, 
        columns=[id_to_db_id[str(i)] for i in range(predictions.shape[1])])

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

    if len(unique_uniprots) > 0:

        # filenames to output enrichment
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

def activity_predict(
    user,
    input_file, 
    compression="zip",
    root_dir="user_files",
    # archive_dir="archives",
    threshold=500,
    enrichment=True):

    '''
    Perform activity prediction with PASS
    '''

    root_dir = os.path.join(root_dir, 
        "user_id={}".format(user.id), "activity_prediction")
    os.makedirs(root_dir, exist_ok=True)

    identifier = determine_identifier(input_file)

    output_dir = os.path.join(root_dir, identifier)
    os.makedirs(output_dir, exist_ok=True)

    input_file = process_input_file(input_file, 
        desired_format=".sdf", output_dir=output_dir)

    base_name, extension = os.path.splitext(input_file)
    assert extension == ".sdf"

    pass_out_file = base_name + "-PASS-out.sdf"

    cmd = "PASS2019toSDF.exe {} {}".format(input_file, pass_out_file)
    print ("executing command:", cmd)

    ret = os.system(cmd)
    assert ret == 0
    assert os.path.exists(pass_out_file)

    # add enrichment to output directory
    if enrichment:
        ret = perform_enrichment_on_PASS_file(pass_out_file,
            output_dir=output_dir, threshold=threshold)
        # delete pass out file
        print ("deleting pass output file:", pass_out_file)
        os.remove(pass_out_file)

    # build zip file containing all targets / run settings / run output
    archive_filename = os.path.join(root_dir,
        identifier)
    print ("writing archive to", 
        archive_filename + "." + compression)

    shutil.make_archive(archive_filename, 
        compression, output_dir)

    attachment_filename = f"{archive_filename}.{compression}"
    send_file_to_user(user, attachment_filename)

    return 0

class User:

    def __init__(self, id):
        self.id = id

if __name__ == "__main__":
    # name = "david"
    # email = "davemcdonald93@gmail.com"
    # input_file = "/home/david/Desktop/test-PASS-out.sdf"

    # # determine_targets(input_file)
    # perform_enrichment_on_PASS_file(input_file, output_dir="/home/david/Desktop", threshold=0)

    # smiles_file = "/home/david/Desktop/pdb_ligands.smi"
    smiles_file = "coconut_smiles.smi"

    # read (and filter smiles)
    smiles = read_smiles(
        smiles_file,
        filter_valid=True, 
        return_series=True)

    n_compounds = smiles.shape[0]

    chunksize = 40000
    n_chunks = n_compounds // chunksize + 1
    for chunk_no in range(n_chunks):
        predictions_filename = f"coconut_uniprot_predictions_chunk_{chunk_no}.csv.gz"
        if os.path.exists(predictions_filename):
            continue
        print ("processing chunk", chunk_no+1)

        chunk = smiles[chunk_no*chunksize:(chunk_no+1)*chunksize]

        predictions = perform_predicton_with_novel_classifier(chunk)

        print ("writing predictions to", predictions_filename)
        predictions.to_csv(predictions_filename)