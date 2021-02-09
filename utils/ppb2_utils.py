
import os
import gzip 
import pickle5 as pkl

import numpy as np
import pandas as pd

from functools import partial

from utils.io import read_smiles, load_json

def load_model(model_filename):
    assert model_filename.endswith(".pkl.gz")
    assert os.path.exists(model_filename), model_filename
    print ("reading model from", model_filename)
    with gzip.open(model_filename, "rb") as f:
        return pkl.load(f)


def perform_predicton_with_novel_classifier(
    smiles,
    # model_filename="models/morg3-xgc.pkl.gz",
    model_filename="models/morg2-nn+nb.pkl.gz",
    k=2000,
    n_proc=1,
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
    model.set_k(k)
    
    # make prediction using pretrained model
    # return as n_targets x n_compounds
    pred = model.predict(smiles).T 
    probs = model.predict_proba(smiles).T 

    # id_to_db_id = load_json("id_to_db_id.json")
    # id_to_target_acc = load_json("models/target_ids.json")
    # id_to_chembl_id = load_json("models/id_to_chembl.json")
    id_to_acc = load_json("models/id_to_uniprot.json")

    pred = pd.DataFrame(pred, 
        index=[id_to_acc[str(i)] 
            for i in range(pred.shape[0])],
        columns=smiles.index, 
    )

    probs = pd.DataFrame(probs, 
        index=[id_to_acc[str(i)] 
            for i in range(probs.shape[0])],
        columns=smiles.index, 
    )

    return pred, probs

def rescale_col(col, max_confidence=1000):
    '''
    linearly scale score by score rank
    '''
    unique_values = sorted(set(col)) # sort unique values in ascending order
    scores = np.linspace(0, max_confidence, len(unique_values)) # divide score evenly by number of unique values
    score_map = {val: score for val, score in zip(unique_values, scores)}
    return col.map(score_map)

def rescale_predicted_uniprot_confidences(predictions, max_confidence=1000):
    assert isinstance(predictions, pd.DataFrame)
    # assert predictions.shape[0] == 1683
    # rescale by max confidence (per compound -- (over all targets)  -- axis 0)
    # return (predictions.divide(predictions.max(axis=0,)) * max_confidence).astype(int)
    # predictions[predictions>0] = max_confidence 
    # return predictions   
    return predictions.apply(partial(rescale_col, max_confidence=max_confidence), axis=0).astype(int)