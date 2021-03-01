import os

import numpy as np
import pandas as pd

from django.shortcuts import render

from django.http import HttpResponse, HttpResponseRedirect
from django.views.static import serve

import html
import urllib.parse as urlparse

from utils.queries import (
    get_all_pathways,
    get_pathway_hits
)

from natural_products.views import (VALID_ORGANISMS, MIN_VAL, MAX_VAL, STEP, 
    DEFAULT_THRESHOLD, MAX_HITS_FOR_IMAGE, filter_columns)

from collections import defaultdict

def all_pathways_view(request):

    organisms = VALID_ORGANISMS

    pathways = get_all_pathways(organisms=organisms)

    pathways = [
        (pathway, urlparse.quote(pathway), organism, urlparse.quote(organism))
            for pathway, organism in pathways
    ]

    context = {
        "pathways": pathways
    }

    return render(request,
        "natural_products/pathways/all_pathways.html",
        context)

def pathway_info_view(request, pathway_organism):

    threshold = DEFAULT_THRESHOLD
    filter_pa_pi = True
    pathway, organism = pathway_organism.split(":_:")
    pathway = urlparse.unquote(pathway)
    organism = urlparse.unquote(organism)

    # query database
    pathways = [pathway]

    pathway_hits, columns = get_pathway_hits(
        pathways, 
        threshold=threshold, filter_pa_pi=filter_pa_pi,
        organism=organism, 
        as_dict=True,
        )
    num_hits = len(pathway_hits)
    show_images = num_hits<MAX_HITS_FOR_IMAGE
    cols_to_remove = {"smiles"}
    if not show_images:
        cols_to_remove.add("image")
    columns = filter_columns(columns, cols_to_remove)

    request.session["targets"] = pathways # for downloading
    request.session["thresholds"] = [threshold]
    request.session["hits"] = pathway_hits
    # request.session["columns"] = columns

    context = {
        "pathway": pathway,
        "organism": organism,
        "threshold": threshold,
        "pathway_hits": pathway_hits,#[:MAX_RECORDS],
        "num_hits": num_hits,
        "columns": columns, 
        "allow_optimise": False
    }

    return render(request,
        "natural_products/pathways/pathway_info.html", 
        context)
        
def pathway_select_view(request):

    organisms = VALID_ORGANISMS

    all_pathways = get_all_pathways(organisms=organisms, )
    pathways = defaultdict(list)
    for p, o in all_pathways:
        pathways[o].append((p, urlparse.quote(p)))

    # pathways = (
    #         (p, urlparse.quote(p), o, urlparse.quote(o))
    #     for p, o in pathways
    # )

    thresholds = range(MAX_VAL, MIN_VAL-1, STEP)
    context = {
        # "organisms": organisms,
        "pathways": dict(pathways),
        "thresholds": thresholds}
    return render(request,
        "natural_products/pathways/pathway_select.html", context)

def make_column_header(s):
    return s.replace("_", " ").title()

def show_pathway_hits_view(request, ):

    organism = request.GET["organism"]
    pathways = request.GET.getlist(f"{organism}-pathways")
    threshold = request.GET["threshold"]
    # filter_pa_pi = request.GET.get("checkbox") == "on"
    filter_pa_pi = True 

    try:
        threshold = int(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")

    # query database
    # organisms = [urlparse.unquote(organism) 
        # for organism in organisms]
    pathways = [urlparse.unquote(pathway) 
        for pathway in pathways]

    pathway_hits, columns = get_pathway_hits(
        pathways, 
        threshold=threshold, 
        filter_pa_pi=filter_pa_pi,
        organism=organism, 
        as_dict=True,
        # limit=100,
        )
    num_hits = len(pathway_hits)
    show_images = num_hits<MAX_HITS_FOR_IMAGE
    cols_to_remove = {"smiles"}
    if not show_images:
        cols_to_remove.add("image")
    columns = filter_columns(columns, cols_to_remove)

    request.session["targets"] = pathways # for downloading
    request.session["threshold"] = threshold
    request.session["hits"] = pathway_hits

    context = {
        "pathways": pathways,
        "threshold": threshold,
        "pathway_hits": pathway_hits,#[:MAX_RECORDS],
        "num_hits": num_hits,
        # "show_images": show_images
        "columns": columns,
        "allow_optimise": False
    }

    return render(request,
        "natural_products/pathways/pathway_hits.html", 
        context)
