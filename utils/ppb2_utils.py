
import os
import gzip 
import pickle5 as pkl

import pandas as pd

from utils.io import read_smiles, load_json

def load_model(model_filename):
    assert model_filename.endswith(".pkl.gz")
    assert os.path.exists(model_filename), model_filename
    print ("reading model from", model_filename)
    with gzip.open(model_filename, "rb") as f:
        return pkl.load(f)


def perform_predicton_with_novel_classifier(
    smiles,
    model_filename="models/morg3-xgc.pkl.gz",
    n_proc=6,
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
    # return as n_targets x n_compounds
    predictions = model.predict_proba(smiles).T 

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
    assert predictions.shape[0] == 1683
    # rescale by max confidence (per compound -- (over all targets)  -- axis 0)
    return (predictions.divide(predictions.max(axis=0,)) * max_confidence).astype(int)