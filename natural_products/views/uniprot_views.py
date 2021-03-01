import os

import numpy as np
import pandas as pd

from django.shortcuts import render

from django.http import HttpResponse, HttpResponseRedirect
from django.views.static import serve

import html
import urllib.parse as urlparse

from utils.queries import (
    get_targets_for_uniprot,
    get_protein_gene_from_acc,
    get_all_uniprot,
    get_all_pathways_for_uniprots,
    get_all_reactions_for_uniprots,
    get_combined_compound_confidences_for_uniprots,
    get_drugs_for_uniprots,
    get_diseases_for_uniprots
)

from natural_products.views import (VALID_ORGANISMS, MIN_VAL, MAX_VAL, STEP, 
    DEFAULT_THRESHOLD, MAX_HITS_FOR_IMAGE, filter_columns)

from collections import defaultdict

def all_accs_view(request):

    # organisms = VALID_ORGANISMS

    uniprot_targets = get_all_uniprot(filter_valid=True)

    context = {
        "uniprot_targets": uniprot_targets
    }

    return render(request,
        "natural_products/uniprot/all_acc_targets.html",
        context)

def acc_info_view(request, acc):

    threshold = DEFAULT_THRESHOLD

    all_info = get_protein_gene_from_acc(acc, as_dict=True)[0]

    # # get inferred NPS
    # inferred_nps, cols = get_inferred_compounds_for_uniprots(acc, threshold=threshold)
    # inferred_nps = [
    #     {c:r}
    #     for record in inferred_nps
    #     for c, r in zip(cols, record)
    # ]

    # # get predicted NPS
    # predicted_nps, cols = get_predicted_compounds_for_uniprots(acc, threshold=threshold)
    # predicted_nps = [
    #     {c:r}
    #     for record in predicted_nps
    #     for c, r in zip(cols, record)
    # ]
    targets, target_cols = get_targets_for_uniprot(acc, as_dict=True)

    nps, np_cols = get_combined_compound_confidences_for_uniprots(acc, threshold=threshold, as_dict=True)
    cols_to_remove = {"smiles"}
    if len(nps) > MAX_HITS_FOR_IMAGE:
        cols_to_remove.add("image")
    np_cols = filter_columns(np_cols, cols_to_remove)

    # get drugs
    drugs, drug_cols = get_drugs_for_uniprots(acc, as_dict=True)

    # get diseases
    diseases, disease_cols = get_diseases_for_uniprots(acc, as_dict=True)

    # # get pathways
    pathways, pathway_cols = get_all_pathways_for_uniprots(acc, as_dict=True )

    # # get reactions
    reactions, reaction_cols = get_all_reactions_for_uniprots(acc, as_dict=True)

    all_data = {
        "AssociatedTargets": (target_cols, targets),
        "NaturalProducts": (np_cols, nps),
        "Drugs": (drug_cols, drugs),
        "Diseases": (disease_cols, diseases),
        "Pathways": (pathway_cols, pathways),
        "Reactions": (reaction_cols, reactions),
    }

    context = {
        "all_info": all_info,
        "threshold": threshold,
        "all_data": all_data,
        # "tables" : {
            # "InferredNPs": inferred_nps,
            # "PredictedNPs": predicted_nps,
        # },
        # "inferred_nps": inferred_nps,#[:MAX_RECORDS],
        # "predicted_nps": predicted_nps,
        # "nps": nps,
        # "drugs": drigs,
        # "diseases": diseases,
        # "pathways": pathways,
        # "reactions": reactions,
        # "show_images": num_hits<MAX_HITS_FOR_IMAGE,
        "allow_optimise": False
    }

    return render(request,
        "natural_products/uniprot/acc_info.html", context)