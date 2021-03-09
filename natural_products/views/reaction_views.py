import os

import numpy as np
import pandas as pd

from django.shortcuts import render

from django.http import HttpResponse, HttpResponseRedirect
from django.views.static import serve

import html
import urllib.parse as urlparse

from utils.queries import (
    get_all_reactions,
    get_reaction_hits
)

from natural_products.views import (VALID_ORGANISMS, THRESHOLDS, DEFAULT_TARGET_COVERAGE,
    DEFAULT_THRESHOLD, MAX_HITS_FOR_IMAGE, filter_columns)

from collections import defaultdict

def all_reactions_view(request):

    organisms = VALID_ORGANISMS

    reactions = get_all_reactions(organisms=organisms)

    reactions = [
        (reaction, urlparse.quote(reaction), organism, urlparse.quote(organism))
            for reaction, organism in reactions
    ]

    context = {
        "reactions": reactions
    }

    return render(request,
        "natural_products/reactions/all_reactions.html",
        context)

def reaction_info_view(request, reaction_organism):

    threshold = DEFAULT_THRESHOLD
    filter_pa_pi = True
    reaction, organism = reaction_organism.split(":_:")
    reaction = urlparse.unquote(reaction)
    organism = urlparse.unquote(organism)

    # query database
    reactions = [reaction]

    reaction_hits, columns = get_reaction_hits(
        reactions, 
        threshold=threshold, filter_pa_pi=filter_pa_pi,
        organism=organism, 
        min_target_coverage=DEFAULT_TARGET_COVERAGE,
        as_dict=True,
        )
    num_hits = len(reaction_hits)
    show_images = num_hits<MAX_HITS_FOR_IMAGE
    cols_to_remove = {"smiles"}
    if not show_images:
        cols_to_remove.add("image")
    columns = filter_columns(columns, cols_to_remove)

    request.session["targets"] = reactions # for downloading
    request.session["thresholds"] = [threshold]
    request.session["hits"] = reaction_hits
    # request.session["columns"] = columns

    context = {
        "reaction": reaction,
        "organism": organism,
        "threshold": threshold,
        "reaction_hits": reaction_hits,#[:MAX_RECORDS],
        "num_hits": num_hits,
        "columns": columns, 
        "allow_optimise": False
    }

    return render(request,
        "natural_products/reactions/reaction_info.html", 
        context)
        
def reaction_select_view(request):

    organisms = VALID_ORGANISMS

    all_reactions = get_all_reactions(organisms=organisms, )
    reactions = dict()
    for p, o in all_reactions:
        if o not in reactions:
            reactions[o] = {"id": o.replace(" ", "_"), "reactions":[]}
        reactions[o]["reactions"].append((p, urlparse.quote(p)))

    context = {
        "reactions": reactions,
        "thresholds": THRESHOLDS
    }
    return render(request,
        "natural_products/reactions/reaction_select.html", context)

def show_reaction_hits_view(request, ):

    organism = request.GET["organism"]
    organism_safe = organism.replace(" ", "_")
    reactions = request.GET.getlist(f"{organism_safe}-reactions")
    threshold = request.GET["threshold"]
    min_reactions_hit = request.GET["min_reactions_hit"]
    min_target_coverage = request.GET["min_target_coverage"]

    # filter_pa_pi = request.GET.get("checkbox") == "on" #TODO
    filter_pa_pi = True

    try:
        threshold = int(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")
    try:
        min_reactions_hits = int(min_reactions_hit)
    except ValueError:
        return HttpResponse("Invalid min_reactions_hit")
    try:
        min_target_coverage = float(min_target_coverage)
    except ValueError:
        return HttpResponse("Invalid min_target_coverage")
    # query database
    # organisms = [urlparse.unquote(organism) 
        # for organism in organisms]
    reactions = [urlparse.unquote(reaction) 
        for reaction in reactions]

    reaction_hits, columns = get_reaction_hits(
        reactions, 
        threshold=threshold, 
        min_target_coverage=min_target_coverage,
        min_reactions_hit=min_reactions_hit,
        filter_pa_pi=filter_pa_pi,
        organism=organism, 
        as_dict=True,
        # limit=100,
        )
    num_hits = len(reaction_hits)
    show_images = num_hits<MAX_HITS_FOR_IMAGE
    cols_to_remove = {"smiles"}
    if not show_images:
        cols_to_remove.add("image")
    columns = filter_columns(columns, cols_to_remove)

    request.session["targets"] = reactions # for downloading
    request.session["threshold"] = threshold
    request.session["hits"] = reaction_hits

    context = {
        "reactions": reactions,
        "threshold": threshold,
        "min_reactions_hit": min_reactions_hit,
        "min_target_coverage": min_target_coverage,
        "reaction_hits": reaction_hits,#[:MAX_RECORDS],
        "num_hits": num_hits,
        # "show_images": show_images
        "columns": columns,
        "allow_optimise": False
    }

    return render(request,
        "natural_products/reactions/reaction_hits.html", 
        context)
