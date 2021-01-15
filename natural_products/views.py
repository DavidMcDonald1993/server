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
    write_smiles_to_file, 
    get_coconut_compound_info_from_mongo,
    draw_molecule, 
)
from utils.queries import (
    get_all_targets_for_categories, 
    get_all_pathways, 
    get_all_reactions, 
    get_all_pathways_for_compounds, 
    get_all_reactions_for_compounds,
    get_target_hits,
    get_pathway_hits,
    get_reaction_hits,
    get_all_activities_for_compound,
    get_info_for_multiple_compounds
)

# Create your views here.

MAX_VAL = 950
MIN_VAL = 650
STEP = -50

MAX_HITS_FOR_IMAGE = 2000

def index_view(request):
    context = {}
    return render(request, 
        "natural_products/index.html", context)
        
def target_select_view(request):

    targets = get_all_targets_for_categories(
        categories={"MECHANISMS", "GENE_EXPRESSION"})
    targets = (
            (c, t, urlparse.quote(t))
        for c, t in targets
    )

    thresholds = range(MAX_VAL, MIN_VAL-1, STEP)

    context = {
        "targets": targets,
        "thresholds": thresholds}
    return render(request,
        "natural_products/target_select.html", context)

def show_target_hits_view(request, ):

    targets = request.GET.getlist("targets")
    threshold = request.GET["threshold"]
    filter_pa_pi = request.GET.get("checkbox") == "on" #TODO

    try:
        threshold = int(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")

    # query database
    targets = [urlparse.unquote(target) 
        for target in targets]

    thresholds = [threshold]
    target_hits, columns = get_target_hits(
        targets, thresholds, filter_pa_pi=filter_pa_pi)
    num_hits = len(target_hits)

    request.session["targets"] = targets
    request.session["thresholds"] = thresholds
    request.session["hits"] = target_hits
    request.session["columns"] = columns

    context = {
        "targets": targets,
        "thresholds": thresholds,
        "target_hits": target_hits,
        "num_hits": num_hits,
        "show_images": num_hits<MAX_HITS_FOR_IMAGE
    }

    return render(request,
        "natural_products/target_hits.html", context)
        
def pathway_select_view(request):

    organisms = ("Homo sapiens", "Rattus norvegicus", "Mus musculus")

    pathways = get_all_pathways(organisms=organisms, )
    pathways = (
            (p, urlparse.quote(p), o, urlparse.quote(o))
        for p, o in pathways
    )

    thresholds = range(MAX_VAL, MIN_VAL-1, STEP)
    context = {
        "pathways": pathways,
        "thresholds": thresholds}
    return render(request,
        "natural_products/pathway_select.html", context)

def show_pathway_hits_view(request, ):

    organisms, pathways = list(zip(*
        (pathway.split(":_:") for pathway in request.GET.getlist("pathways"))))
    threshold = request.GET["threshold"]
    filter_pa_pi = request.GET.get("checkbox") == "on" 
    # organism = "Homo sapiens"

    try:
        threshold = int(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")

    # query database
    organisms = [urlparse.unquote(organism) 
        for organism in organisms]
    pathways = [urlparse.unquote(pathway) 
        for pathway in pathways]

    pathway_hits, columns = get_pathway_hits(
        pathways, 
        threshold=threshold, filter_pa_pi=filter_pa_pi,
        organisms=organisms, 
        # limit=100,
        )
    num_hits = len(pathway_hits)

    request.session["targets"] = pathways # for downloading
    request.session["thresholds"] = [threshold]
    request.session["hits"] = pathway_hits
    request.session["columns"] = columns

    # information requested for each pathway 
    pathway_columns = (
        "targets", "num_targets", 
        "accs", "num_accs", 
        "total_accs", "coverage", 
        "pathway_name", "organism", "url"
    )
    num_cols = len(pathway_columns)
    num_pathways = len(pathways)

    pathway_hits = [
        
        {
            "id": compound_id, 
            "image": image,
            "name": compound_name, 
            "formula": compound_formula, 
            "smiles": compound_smiles,
            "pathways": [
                {   
                    col: pathway_values[pathway_number*num_cols+col_num]
                        for col_num, col in enumerate(pathway_columns)} 
                for pathway_number in range(num_pathways)   
            ]
        }
        for (compound_id, image, 
            compound_name, compound_formula, compound_smiles,
            *pathway_values) in pathway_hits
    ]

    # from decimal import Decimal

    # for hit in pathway_hits[:1]:
    #     for key in hit:
    #         # assert not isinstance(hit[key], Decimal)
    #         print (key, hit[key], type(hit[key]))

    #     for pathway in hit["pathways"]:
    #         for key in pathway:
    #             # assert not isinstance(pathway[key], Decimal)
    #             print (key, pathway[key], type(pathway[key]))

    context = {
        "pathways": pathways,
        "threshold": threshold,
        "pathway_hits": pathway_hits,
        "num_hits": num_hits,
        "show_images": num_hits<MAX_HITS_FOR_IMAGE
    }

    return render(request,
        "natural_products/pathway_hits.html", 
        context)

def reaction_select_view(request):

    organisms = ("Homo sapiens", "Rattus norvegicus", "Mus musculus")

    reactions = get_all_reactions(organisms=organisms)
    reactions = (
            (r, urlparse.quote(r), o, urlparse.quote(o))
        for r, o in reactions
    )

    thresholds = range(MAX_VAL, MIN_VAL-1, STEP)
    context = {
        "reactions": reactions,
        "thresholds": thresholds,
    }
    return render(request,
        "natural_products/reaction_select.html", 
        context)

def show_reaction_hits_view(request, ):

    organisms, reactions = list(zip(*
        (reaction.split(":_:") for reaction in request.GET.getlist("reactions"))))
    threshold = request.GET["threshold"]
    filter_pa_pi = request.GET.get("checkbox") == "on"

    try:
        threshold = int(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")

    # query database
    organisms = [urlparse.unquote(organism)
        for organism in organisms]
    reactions = [urlparse.unquote(reaction) 
        for reaction in reactions]

    # thresholds = [threshold]
    reaction_hits, columns = get_reaction_hits(
        reactions, 
        threshold=threshold, filter_pa_pi=filter_pa_pi,
        organisms=organisms,
        # limit=100
        )
    num_hits = len(reaction_hits)

    request.session["targets"] = reactions # for downloading
    request.session["thresholds"] = [threshold]
    request.session["hits"] = reaction_hits
    request.session["columns"] = columns

        # information requested for each pathway 
    reaction_columns = (
        "targets", "num_targets", 
        "accs", "num_accs", 
        "total_accs", "coverage",
        "reaction_name", "organism", "url"
    )
    num_cols = len(reaction_columns)
    num_reactions = len(reactions)

    reaction_hits = [
        {
            "id": compound_id, 
            "image": image,
            "name": compound_name, 
            "formula": compound_formula, 
            "smiles": compound_smiles,
            "reactions": [
                {   
                    col: reaction_values[reaction_number*num_cols+col_num]
                        for col_num, col in enumerate(reaction_columns)} 
                for reaction_number in range(num_reactions)   
            ]
        }
        for compound_id, image, compound_name, compound_formula, compound_smiles,
            *reaction_values in reaction_hits
    ]

    context = {
        "reactions": reactions,
        "threshold": threshold,
        "reaction_hits": reaction_hits,
        "num_hits": num_hits,
        "show_images": num_hits<MAX_HITS_FOR_IMAGE
    }

    return render(request,
        "natural_products/reaction_hits.html", context)


def all_compounds_view(request):

    compounds = get_info_for_multiple_compounds()
    context = {
        "compounds": compounds
    }

    return render(request, 
        "natural_products/all_compounds.html", context)

def compound_info_view(request, compound_id):

    compound_id = "CNP" + compound_id

    threshold = 750

    # query database
    compound_info = get_coconut_compound_info_from_mongo(compound_id)

    context = {
        "threshold": threshold
    }

    assert "clean_smiles" in compound_info
    assert "name" in compound_info 
    compound_name = compound_info["name"]

    img_filename = get_info_for_multiple_compounds(
        compound_id, columns=("image_path", ))[0][0]
    if os.path.exists(os.path.join("static", img_filename)):
        context["img_filename"] = img_filename

    activities = get_all_activities_for_compound(
        compound_id, 
        threshold=threshold, 
        )


    pathways = get_all_pathways_for_compounds(
        compound_id, 
        threshold=threshold,
        # organism="Homo sapiens"
        )

    reactions = get_all_reactions_for_compounds(
        compound_id,  
        threshold=threshold,
        # organism="Homo sapiens"
        )

    context.update({
        "compound_id": compound_id,
        "compound_name": compound_name,
        "compound_info": compound_info,
        "activities": activities,
        "pathways": pathways,
        "reactions": reactions,
    })

    return render(request,
        "natural_products/compound_info.html", 
        context)

def download_hits_view(request):

    assert request.user.is_authenticated
    
    user_id = request.user.id
    targets = request.session["targets"]
    thresholds = request.session["thresholds"]
    hits = request.session["hits"]
    columns = request.session["columns"]

    hits = pd.DataFrame(hits, columns=columns)

    record_filename = write_records_to_file(user_id, targets, thresholds, hits)

    return serve(request, 
            os.path.basename(record_filename), 
            os.path.dirname(record_filename))

def optimise_target_hits_view(request):

    assert request.user.is_authenticated
    user_id = request.user.id
    targets = request.session["targets"]
    thresholds = request.session["thresholds"]
    hits = request.session["hits"]

    smiles = get_multiple_compound_info(
        compounds=[hit[0] for hit in hits],
        columns=("coconut_id", "clean_smiles"))

    smiles_filename = write_smiles_to_file(user_id, targets, thresholds, smiles)
    request.session["smiles_filename"] = smiles_filename

    return HttpResponseRedirect("/hit_optimisation/upload")