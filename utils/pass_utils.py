import sys
import os.path
sys.path.insert(1, 
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from utils.io import load_json
from utils.rdkit_utils import LoadSDF

import pandas as pd

from collections import defaultdict

def get_categories():
    return sorted({
        "EFFECTS",
        "MECHANISMS",
        'TOXICITY', 
        'ANTITARGETS',
        'METABOLISM', 
        'GENE_EXPRESSION', 
        'TRANSPORTERS'
    })

def get_all_compounds():
    return load_json("compound_ids.json")

def get_all_targets():
    return load_json("target_ids.json")

def get_targets_for_category(category):
    assert category in get_categories()
    return load_json(f"{category}_targets.json",)

def remove_invalid_characters(s):
    return s.replace(".", "") # also allows use of INT storage in database

def parse_pass_spectra(s, mapping=None):
    split = s.split()
    pa = int(split[0])
    pi = int(split[1])
    target = " ".join(split[2:])
    if mapping is not None:
        assert target in mapping 
        target = mapping[target] # map to ids
    return target, {"Pa": pa, "Pi": pi, "Conf": pa-pi}

def parse_pass_confidence(s, mapping=None):
    split = s.split()
    pa = int(split[0])
    pi = int(split[1])
    target = " ".join(split[2:])
    if mapping is not None:
        assert target in mapping 
        target = mapping[target] # map to ids
    return target, pa-pi

def determine_targets(prediction_SDF_filename, threshold=0):
    '''
    determine targets from PASS-predicted SDF file (for enrichment analysis)
    only valid for cmd version
    '''
    print ("reading PASS predictions from", prediction_SDF_filename)
    assert os.path.exists(prediction_SDF_filename) 

    prediction_chunks = LoadSDF(prediction_SDF_filename)

    active_targets = defaultdict(list)

    for chunk in prediction_chunks:
        if "PASS_ACTIVITY_SPECTRUM" not in chunk.columns:
            continue
        if not (chunk["ID"]=="").any():
            chunk = chunk.set_index("ID", drop=True)
        chunk = chunk.loc[~pd.isnull(chunk["PASS_ACTIVITY_SPECTRUM"])]

        chunk_activities = chunk["PASS_ACTIVITY_SPECTRUM"].map(
            lambda s:
            map(parse_pass_spectra, 
                map(remove_invalid_characters, s.split("\n"))),
            na_action="ignore")

        for compound, target_activities in chunk_activities.items():
            for target, activities in target_activities:
                if activities["Pa"] > threshold and activities["Pa"] > activities["Pi"]:
                    active_targets[compound].append((target, activities["Conf"]))

    return active_targets
    
def determine_confidences(prediction_SDF_filename, threshold=0):
    '''
    determine confidences from PASS-predicted SDF file (for enrichment analysis)
    only valid for cmd version
    '''
    print ("reading PASS predictions from", prediction_SDF_filename)
    assert os.path.exists(prediction_SDF_filename) 

    prediction_chunks = LoadSDF(prediction_SDF_filename)

    all_confidences = dict()#pd.DataFrame(dtype=int)

    for chunk in prediction_chunks:
        if "PASS_ACTIVITY_SPECTRUM" not in chunk.columns:
            continue
        if not (chunk["ID"]=="").any():
            chunk = chunk.set_index("ID", drop=True)
        chunk = chunk.loc[~pd.isnull(chunk["PASS_ACTIVITY_SPECTRUM"])]

        chunk_confidences = chunk["PASS_ACTIVITY_SPECTRUM"].map(
            lambda s:
                map(parse_pass_confidence, 
                    map(remove_invalid_characters, s.split("\n"))),
                na_action="ignore")

        for compound, target_confidences in chunk_confidences.items():
            targets, confidences = zip(*target_confidences)
            all_confidences[compound] = pd.Series(confidences, index=targets)
            # for target, confidence in confidences:
            #     if activities["Pa"] > threshold and activities["Pa"] > activities["Pi"]:
            #         active_targets[compound].append((target, activities["Conf"]))

    return pd.DataFrame(all_confidences)
        
if __name__ == "__main__":
    # for category in get_categories():
    #     get_targets_for_category(category)

    # get_all_targets()

    name = "david"
    email = "davemcdonald93@gmail.com"
    input_file = "/home/david/Desktop/test-PASS-out.sdf"

    # active_targets = determine_targets(input_file, threshold=500)
    # print (active_targets)
    confidences = determine_confidences(input_file)

    for threshold in [500]:

        # thresholded_confidences = confidences > threshold
        thresholded_confidences = {compound: 
            [k for k, v in confidences[compound].items() if v>threshold] 
            for compound in confidences}


        print (thresholded_confidences)

    # max_confidences = confidences.max(axis=1)

    # # print(max_confidences)

    # from utils.queries import get_uniprots_for_targets

    # records = get_uniprots_for_targets(max_confidences.index)

    # df = pd.DataFrame([ (target, max_confidences[target], uniprot, association_score)
    #     for target, uniprot, association_score in records],
    #     columns=["target", "max_confidence", "ACC", "association_score"])

    # print (df.head())

