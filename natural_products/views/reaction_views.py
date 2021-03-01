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

from natural_products.views import (VALID_ORGANISMS, MIN_VAL, MAX_VAL, STEP, 
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
    reactions = defaultdict(list)
    for p, o in all_reactions:
        reactions[o].append((p, urlparse.quote(p)))

    # reactions = (
    #         (p, urlparse.quote(p), o, urlparse.quote(o))
    #     for p, o in reactions
    # )

    thresholds = range(MAX_VAL, MIN_VAL-1, STEP)
    context = {
        # "organisms": organisms,
        "reactions": dict(reactions),
        "thresholds": thresholds}
    return render(request,
        "natural_products/reactions/reaction_select.html", context)

def make_column_header(s):
    return s.replace("_", " ").title()

def show_reaction_hits_view(request, ):

    organism = request.GET["organism"]
    reactions = request.GET.getlist(f"{organism}-reactions")
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
    reactions = [urlparse.unquote(reaction) 
        for reaction in reactions]

    reaction_hits, columns = get_reaction_hits(
        reactions, 
        threshold=threshold, 
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
        "reaction_hits": reaction_hits,#[:MAX_RECORDS],
        "num_hits": num_hits,
        # "show_images": show_images
        "columns": columns,
        "allow_optimise": False
    }

    return render(request,
        "natural_products/reactions/reaction_hits.html", 
        context)
