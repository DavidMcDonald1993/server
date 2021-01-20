
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
    get_all_targets, determine_targets, determine_confidences)
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
    # max_confidence=1000
    ):
    '''
    Predict from SMILES using novel classifier
    '''
    if isinstance(smiles, str): # read from file
        assert smiles.endswith(".smi")
    
        # read (and filter smiles)
        smiles = read_smiles(
            smiles,
            filter_valid=True, 
            return_series=True)

    assert isinstance(smiles, pd.Series)

    model = load_model(model_filename)
    if hasattr(model, "n_proc"):
        model.set_n_proc(n_proc)
    
    # make prediction using pretrained model
    predictions = model.predict_proba(smiles).T 

    # return as n_targets x n_compounds

    # id_to_db_id = load_json("id_to_db_id.json")
    id_to_target_acc = load_json("models/target_ids.json")
    
    return pd.DataFrame(predictions, 
        # index=[id_to_db_id[str(i)] for i in range(predictions.shape[1])],
        index=[id_to_target_acc[str(i)] 
            for i in range(predictions.shape[0])],
        columns=smiles.index, 
    )

def rescale_predicted_uniprot_confidences(predictions, max_confidence=1000):
    assert isinstance(predictions, pd.DataFrame)
    # rescale by max confidence (per compound -- axis 1)
    return (predictions.divide(predictions.max(axis=1,), axis=0) * max_confidence).astype(int)

def perform_enrichment_on_uniprot_accs(
    targets_to_uniprot, 
    output_dir,
    threshold=500,
    ):

    output_dir = os.path.join(output_dir, "enrichment")
    os.makedirs(output_dir, exist_ok=True)

    print ("perfoming enrichment analysis on targets_to_uniprot file",
        "to directory", output_dir)

    # # determine active targets from PASS-predicted SDF file
    # active_targets = determine_targets(pass_out_file, threshold=threshold)
    # active_targets_filename = os.path.join(output_dir,
    #     f"active_targets_threshold={threshold}.json")
    # write_json(active_targets, active_targets_filename)

    # unique_target_names = {target 
    #     for compound, targets in active_targets.items()
    #     for target, confidence in targets}
    
    # # unique_uniprots = targets_to_uniprot_ids(unique_target_names)
    # targets_to_uniprot = get_uniprots_for_targets(unique_target_names) # SQL query
    # targets_to_uniprot_filename = os.path.join(output_dir, "targets_to_uniprot.csv")
    # print ("writing targets to uniprot to", targets_to_uniprot_filename)
    # targets_to_uniprot = pd.DataFrame(targets_to_uniprot, 
    #     columns=["target", "uniprot_ACC", "association_score"])
    # targets_to_uniprot.to_csv(targets_to_uniprot_filename)

    above_threshold = targets_to_uniprot.loc[
        targets_to_uniprot["max_confidence"] > threshold]

    unique_uniprots = set(above_threshold["uniprot_ACC"])
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

def write_actives(confidence_df, threshold, output_dir):
    actives = {compound:
        [target for target, confidence in confidence_df[compound].items()
            if confidence > threshold]
            for compound in confidence_df}
    active_targets_filename = os.path.join(output_dir,
        f"actives_(threshold={threshold}).json")
    write_json(actives, active_targets_filename)

def activity_predict(
    user,
    input_file, 
    pass_predict=True,
    ppb2_predict=True,
    compression="zip",
    root_dir="user_files",
    # archive_dir="archives",
    threshold=500,
    enrichment=True):

    '''
    Perform activity prediction with PASS and PPB2
    '''

    root_dir = os.path.join(root_dir, 
        "user_id={}".format(user.id), "activity_prediction")
    os.makedirs(root_dir, exist_ok=True)

    identifier = determine_identifier(input_file)

    output_dir = os.path.join(root_dir, identifier)
    os.makedirs(output_dir, exist_ok=True)

    if pass_predict:

        '''
        BEGIN PASS
        '''

        pass_output_dir = os.path.join(output_dir, "PASS")
        os.makedirs(pass_output_dir, exist_ok=True)

        input_sdf_file = process_input_file(input_file, 
            desired_format=".sdf", output_dir=pass_output_dir)

        base_name, extension = os.path.splitext(input_sdf_file)
        assert extension == ".sdf"

        pass_out_file = base_name + "-PASS-out.sdf"

        cmd = "PASS2019toSDF.exe {} {}".format(input_sdf_file, pass_out_file)
        print ("executing command:", cmd)

        ret = os.system(cmd)
        assert ret == 0
        assert os.path.exists(pass_out_file)

        # DF n_targets * n_compounds
        confidences = determine_confidences(pass_out_file)

        # delete large SDF file
        print ("deleting pass output file:", pass_out_file)
        os.remove(pass_out_file)

        # write confidences to file
        confidence_filename = os.path.join(pass_output_dir, 
            "unthresholded_confidences.csv")
        print ("writing confidences to", confidence_filename)
        confidences.to_csv(confidence_filename)

        # write jsons for threshold(s)?
        for threshold in range(500, 1000, 100):
            # active_targets = {compound:
            #     [target for target, confidence in confidences[compound].items()
            #         if confidence > threshold]
            #         for compound in confidences}
            # active_targets_filename = os.path.join(pass_output_dir,
            #     f"active_targets_threshold={threshold}.json")
            # write_json(active_targets, active_targets_filename)
            write_actives(confidences, threshold, output_dir=pass_output_dir)

        max_confidences = confidences.max(axis=1)
        # threshold?

        # get uniprots for targets
        targets_to_uniprot = get_uniprots_for_targets(max_confidences.index)

        targets_to_uniprot_df = pd.DataFrame(
            sorted([ (target, max_confidences[target], uniprot, association_score)
            for target, uniprot, association_score in targets_to_uniprot], 
                key=lambda x: x[1], reverse=True), # sort by confidence
            columns=["target", "max_confidence", "uniprot_ACC", "association_score"])
        targets_to_uniprot_filename = os.path.join(pass_output_dir,
            "unthresholded_targets_to_uniprot.csv")
        targets_to_uniprot_df.to_csv(targets_to_uniprot_filename)
    
    if ppb2_predict:

        '''
        BEGIN PPB2
        '''

        ppb2_output_dir = os.path.join(output_dir, "PPB2")
        os.makedirs(ppb2_output_dir, exist_ok=True)

        input_smiles_file = process_input_file(input_file, 
            desired_format=".smi", output_dir=ppb2_output_dir)

        # use classifier to generate uniprot confidences
        novel_classifier_predictions = perform_predicton_with_novel_classifier(
            input_smiles_file)

        novel_classifier_prediction_filename = os.path.join(ppb2_output_dir,
            "PPB2_uniprot_predictions.csv")
        print ("writing novel classifier predictions to", 
            novel_classifier_prediction_filename)
        novel_classifier_predictions.to_csv(novel_classifier_prediction_filename)

        # rescale by compound to range [0, 1000]
        novel_classifier_predictions = \
            rescale_predicted_uniprot_confidences(novel_classifier_predictions)

        novel_classifier_prediction_filename = os.path.join(ppb2_output_dir,
            "PPB2_uniprot_predictions_rescaled.csv")
        print ("writing novel classifier predictions to", 
            novel_classifier_prediction_filename)
        novel_classifier_predictions.to_csv(novel_classifier_prediction_filename)

        # write jsons for threshold(s)?
        for threshold in range(500, 1000, 100):
            # active_targets = {compound:
            #     [target for target, confidence in novel_classifier_predictions[compound].items()
            #         if confidence > threshold]
            #         for compound in novel_classifier_predictions}
            # active_targets_filename = os.path.join(ppb2_output_dir,
            #     f"active_uniprots_threshold={threshold}.json")
            # write_json(active_targets, active_targets_filename)
            write_actives(novel_classifier_predictions, threshold, 
                output_dir=ppb2_output_dir)

    if enrichment:

        '''
        BEGIN ENRICHMENT ANALYSIS TODO incorporate PPB2?
        '''

        ret = perform_enrichment_on_uniprot_accs(
            targets_to_uniprot_df,
            output_dir=output_dir, threshold=threshold)
    

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

    def __init__(self, id_):
        self.id = id_

if __name__ == "__main__":
    # name = "david"
    # email = "davemcdonald93@gmail.com"
    # input_file = "/home/david/Desktop/test-PASS-out.sdf"
    user = User(1)
    input_file = "/home/david/Desktop/test.txt"

    ret = activity_predict(user,
        input_file,
        pass_predict=False,
        enrichment=False)

    # # determine_targets(input_file)
    # perform_enrichment_on_PASS_file(input_file, output_dir="/home/david/Desktop", threshold=0)

    # smiles_file = "/home/david/Desktop/pdb_ligands.smi"
    # smiles_file = "coconut_smiles.smi"

    # # read (and filter smiles)
    # smiles = read_smiles(
    #     smiles_file,
    #     filter_valid=True, 
    #     return_series=True)

    # n_compounds = smiles.shape[0]

    # chunksize = 40000
    # n_chunks = n_compounds // chunksize + 1
    # for chunk_no in range(n_chunks):
    #     predictions_filename = f"coconut_uniprot_predictions_chunk_{chunk_no}.csv.gz"
    #     if os.path.exists(predictions_filename):
    #         continue
    #     print ("processing chunk", chunk_no+1)

    #     chunk = smiles[chunk_no*chunksize:(chunk_no+1)*chunksize]

    #     predictions = perform_predicton_with_novel_classifier(chunk)

    #     print ("writing predictions to", predictions_filename)
    #     predictions.to_csv(predictions_filename)