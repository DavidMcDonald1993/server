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

from natural_products.views import (VALID_ORGANISMS, THRESHOLDS, DEFAULT_TARGET_COVERAGE,
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
        min_target_coverage=DEFAULT_TARGET_COVERAGE,
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
    pathways = dict()
    for p, o in all_pathways:
        if o not in pathways:
            pathways[o] = {"id": o.replace(" ", "_"), "pathways":[]}
        pathways[o]["pathways"].append((p, urlparse.quote(p)))

    context = {
        "pathways": pathways,
        "thresholds": THRESHOLDS
    }
    return render(request,
        "natural_products/pathways/pathway_select.html", context)

def show_pathway_hits_view(request, ):

    organism = request.GET["organism"]
    organism_safe = organism.replace(" ", "_")
    pathways = request.GET.getlist(f"{organism_safe}-pathways")
    threshold = request.GET["threshold"]
    min_pathways_hit = request.GET["min_pathways_hit"]
    min_target_coverage = request.GET["min_target_coverage"]

    # filter_pa_pi = request.GET.get("checkbox") == "on" #TODO
    filter_pa_pi = True

    try:
        threshold = int(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")
    try:
        min_pathways_hits = int(min_pathways_hit)
    except ValueError:
        return HttpResponse("Invalid min_pathways_hit")
    try:
        min_target_coverage = float(min_target_coverage)
    except ValueError:
        return HttpResponse("Invalid min_target_coverage")
    # query database
    # organisms = [urlparse.unquote(organism) 
        # for organism in organisms]
    pathways = [urlparse.unquote(pathway) 
        for pathway in pathways]

    pathway_hits, columns = get_pathway_hits(
        pathways, 
        threshold=threshold, 
        min_target_coverage=min_target_coverage,
        min_pathways_hit=min_pathways_hit,
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
        "min_pathways_hit": min_pathways_hit,
        "min_target_coverage": min_target_coverage,
        "pathway_hits": pathway_hits,#[:MAX_RECORDS],
        "num_hits": num_hits,
        # "show_images": show_images
        "columns": columns,
        "allow_optimise": False
    }

    return render(request,
        "natural_products/pathways/pathway_hits.html", 
        context)
