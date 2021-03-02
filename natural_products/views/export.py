import os

import numpy as np
import pandas as pd

from django.shortcuts import render

from django.http import HttpResponse, HttpResponseRedirect
from django.views.static import serve

import html
import urllib.parse as urlparse

from natural_products.backend import (
    write_records_to_file, 
    write_smiles_to_file
)

def export_hits_view(request):

    if not request.user.is_authenticated:
        return HttpResponseRedirect("/login")
    
    user_id = request.user.id
    targets = request.session["targets"]
    threshold = request.session["threshold"]
    hits = request.session["hits"]
    assert isinstance(hits, list)
    # columns = request.session["columns"]

    hits = pd.DataFrame(hits, )

    record_filename = write_records_to_file(user_id, targets, threshold, hits)

    return serve(request, 
            os.path.basename(record_filename), 
            os.path.dirname(record_filename))

def export_hit_smiles_view(request):

    if not request.user.is_authenticated:
        return HttpResponseRedirect("/login")

    user_id = request.user.id
    targets = request.session["targets"]
    threshold = request.session["threshold"]
    hits = request.session["hits"]
    assert isinstance(hits, list)
    # columns = request.session["columns"]

    # assert "SMILES" in columns
    # assert "ID" in columns 

    # hits = pd.DataFrame(hits, columns=columns)

    smiles = [(hit["id"], hit["smiles"])
        for hit in hits]

    smiles_filename = write_smiles_to_file(user_id, targets, threshold, smiles)

    return serve(request, 
            os.path.basename(smiles_filename), 
            os.path.dirname(smiles_filename))

def optimise_target_hits_view(request):

    if not request.user.is_authenticated:
        return HttpResponseRedirect("/login")
    user_id = request.user.id
    if "targets" not in request.session.keys():
        return HttpResponseRedirect("/")
        
    targets = request.session["targets"]
    threshold = request.session["threshold"]
    hits = request.session["hits"]
    assert isinstance(hits, list)
    # columns = request.session["columns"]

    # assert "SMILES" in columns
    # assert "ID" in columns 

    # hits = pd.DataFrame(hits, columns=columns)

    smiles = [(hit["id"], hit["smiles"])
        for hit in hits]

    smiles_filename = write_smiles_to_file(user_id, targets, threshold, smiles)
    request.session["smiles_filename"] = smiles_filename

    request.session["optimise"] = True

    return HttpResponseRedirect("/hit_optimisation")