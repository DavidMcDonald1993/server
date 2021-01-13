
import os
import gzip 
import pickle5 as pkl

def load_model(model_filename):
    assert model_filename.endswith(".pkl.gz")
    assert os.path.exists(model_filename), model_filename
    print ("reading model from", model_filename)
    with gzip.open(model_filename, "rb") as f:
        return pkl.load(f)