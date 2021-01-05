import os

import numpy as np

from django.shortcuts import render

from django.http import HttpResponse, HttpResponseRedirect
from django.views.static import serve

import html
import urllib.parse as urlparse

from natural_products.backend import (write_records_to_file, write_smiles_to_file, 
    query_target_hits, get_compound_info, draw_molecule, get_multiple_compound_info,
    query_pathway_hits, query_reaction_hits)

from utils.mysql_utils import (get_all_targets_and_categories, get_all_pathways, get_all_reactions)

# Create your views here.

# with open("PASS_targets.txt", "r") as f:
#     pass_targets = [line.rstrip()
#         for line in f.readlines()]

def index(request):
    context = {}
    return render(request, 
        "natural_products/index.html", context)
        
def target_select(request):

    targets = get_all_targets_and_categories()
    targets = (
            (c, t, urlparse.quote(t))
        for c, t in targets
    )

    thresholds = ["{:.02f}".format(threshold)
        for threshold in np.arange(0, 1, 0.05)[::-1]]
    context = {
        "targets": targets,
        "thresholds": thresholds}
    return render(request,
        "natural_products/target_select.html", context)

def show_target_hits(request, ):

    targets = request.GET.getlist("targets")
    threshold = request.GET["threshold"]
    filter_pa_pi = request.GET.get("checkbox") == "on" #TODO
    # filter_pa_pi = False

    try:
        threshold = int(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")

    # query database
    targets = [urlparse.unquote(target) 
        for target in targets]

    thresholds = [threshold]
    target_hits = query_target_hits(targets, thresholds, filter_pa_pi=filter_pa_pi)

    request.session["targets"] = targets
    request.session["thresholds"] = thresholds
    request.session["target_hits"] = target_hits

    context = {
        "targets": targets,
        "thresholds": thresholds,
        "target_hits": target_hits,
        "num_hits": len(target_hits)
    }

    return render(request,
        "natural_products/target_hits.html", context)
        
def pathway_select(request):

    organism = "Homo sapiens"

    pathways = get_all_pathways(organism=organism)
    pathways = (
            (p, urlparse.quote(p), o, urlparse.quote(o))
        for p, o in pathways
    )

    thresholds = ["{:.02f}".format(threshold)
        for threshold in np.arange(0, 1, 0.05)[::-1]]
    context = {
        "pathways": pathways,
        "thresholds": thresholds}
    return render(request,
        "natural_products/pathway_select.html", context)

def show_pathway_hits(request, ):

    pathways = request.GET.getlist("pathways")
    threshold = request.GET["threshold"]
    filter_pa_pi = request.GET.get("checkbox") == "on" #TODO
    # filter_pa_pi = False
    organism = "Homo sapiens"

    try:
        threshold = int(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")

    # query database
    pathways = [urlparse.unquote(pathway) 
        for pathway in pathways]

    # thresholds = [threshold]
    pathway_hits = query_pathway_hits(pathways, 
        threshold=threshold, filter_pa_pi=filter_pa_pi,
        organism=organism, 
        limit=100,
        )

    request.session["pathways"] = pathways
    request.session["threshold"] = threshold
    request.session["pathway_hits"] = pathway_hits

    context = {
        "pathways": pathways,
        "threshold": threshold,
        "pathway_hits": pathway_hits,
        "num_hits": len(pathway_hits)
    }

    return render(request,
        "natural_products/pathway_hits.html", context)

def reaction_select(request):

    organism = "Homo sapiens"

    reactions = get_all_reactions(organism=organism)
    reactions = (
            (r, urlparse.quote(r), o, urlparse.quote(o))
        for r, o in reactions
    )

    thresholds = ["{:.02f}".format(threshold)
        for threshold in np.arange(0, 1, 0.05)[::-1]]
    context = {
        "reactions": reactions,
        "thresholds": thresholds}
    return render(request,
        "natural_products/reaction_select.html", context)

def show_reaction_hits(request, ):

    reactions = request.GET.getlist("reactions")
    threshold = request.GET["threshold"]
    filter_pa_pi = request.GET.get("checkbox") == "on" #TODO
    # filter_pa_pi = False
    organism = "Homo sapiens"

    try:
        threshold = int(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")

    # query database
    reactions = [urlparse.unquote(reaction) 
        for reaction in reactions]

    # thresholds = [threshold]
    reaction_hits = query_reaction_hits(reactions, 
        threshold=threshold, filter_pa_pi=filter_pa_pi,
        organism=organism,
        limit=100)

    request.session["reactions"] = reactions
    request.session["threshold"] = threshold
    request.session["reaction_hits"] = reaction_hits

    context = {
        "reactions": reactions,
        "threshold": threshold,
        "reaction_hits": reaction_hits,
        "num_hits": len(reaction_hits)
    }

    return render(request,
        "natural_products/reaction_hits.html", context)

# def pathway_enrichment(request):

#     context = {}

#     return render(request,
#         "natural_products/pathway_enrichment.html", context)

def all_compounds(request):

    # get compound data from database
    compounds = get_multiple_compound_info() # returns list of tuples
    context = {"compounds": compounds}

    return render(request, 
        "natural_products/all_compounds.html", context)

def compound_info(request, compound_id):

    compound_id = "CNP" + compound_id

    # query database
    info, activities = get_compound_info(compound_id)

    context = {}

    assert isinstance(info, dict), type(info)
    assert "clean_smiles" in info
    assert "name" in info 
    compound_name = info["name"]

    smiles = info["clean_smiles"]
    if smiles is not None:
        img_filename = draw_molecule(smiles)
        if img_filename is not None:
            context["img_filename"] = img_filename

    # convert to list of dicts for easy presentation
    info = info.items()

    context.update({
        "compound_id": compound_id,
        "compound_name": compound_name,
        "info": info,
        "activities": activities})

    return render(request,
        "natural_products/compound_info.html", context)

def download_target_hits(request):

    assert request.user.is_authenticated
    
    username = request.user.username
    targets = request.session["targets"]
    thresholds = request.session["thresholds"]
    records = request.session["records"]

    record_filename = write_records_to_file(username, targets, thresholds, records)

    return serve(request, 
            os.path.basename(record_filename), 
            os.path.dirname(record_filename))

def optimise_target_hits(request):

    assert request.user.is_authenticated
    # assert request.session["category"] == "GENE_EXPRESSION"
    username = request.user.username
    records = request.session["records"]

    smiles = get_multiple_compound_info(
        compounds=(record[0] for record in records),
        columns=("coconut_id", "clean_smiles"))

    smiles_filename = write_smiles_to_file(username, smiles)
    request.session["smiles_filename"] = smiles_filename

    return HttpResponseRedirect("/hit_optimisation/upload")