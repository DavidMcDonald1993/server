
import os

import numpy as np
import pandas as pd

from django.shortcuts import render

from django.http import HttpResponse, HttpResponseRedirect
from django.views.static import serve

import html
import urllib.parse as urlparse

from natural_products.views import MIN_VAL, MAX_VAL, MAX_HITS_FOR_IMAGE, STEP, DEFAULT_THRESHOLD, filter_columns

from utils.queries import get_all_targets_for_categories, get_target_hits

def all_targets_view(request):

    targets, columns = get_all_targets_for_categories(as_dict=False)

    targets = (
            (c, t, urlparse.quote(t))
        for c, t in targets
    )

    context = {
        "targets": targets,
        "columns": columns
    }

    return render(request, 
        "natural_products/targets/all_targets.html",
        context)

def target_info_view(request, target):

    threshold = DEFAULT_THRESHOLD
    filter_pa_pi = True
    target = urlparse.unquote(target)

    # query database
    targets = [target]

    target_hits, columns = get_target_hits(
        targets, threshold, filter_pa_pi=filter_pa_pi,
        as_dict=True)
    num_hits = len(target_hits)
    show_images = num_hits<MAX_HITS_FOR_IMAGE
    cols_to_remove = {"smiles"}
    if not show_images:
        cols_to_remove.add("image")
    columns = filter_columns(columns, cols_to_remove)

    request.session["targets"] = targets
    request.session["threshold"] = threshold
    request.session["hits"] = target_hits


    context = {
        "target": target,
        "threshold": threshold,
        "target_hits": target_hits,
        "num_hits": num_hits,
        "columns": columns,
        "allow_optimise": True
    }

    return render(request, 
        "natural_products/targets/target_info.html",
        context)

def target_select_view(request):

    targets, cols = get_all_targets_for_categories(
        as_dict=False,
        categories={"MECHANISMS", "GENE_EXPRESSION", "TOXICITY"})
    targets = (
            (c, t, urlparse.quote(t))
        for c, t in targets
    )

    thresholds = range(MAX_VAL, MIN_VAL-1, STEP)

    context = {
        "targets": targets,
        "thresholds": thresholds}
    return render(request,
        "natural_products/targets/target_select.html", 
        context)

def show_target_hits_view(request, ):

    targets = request.GET.getlist("targets")
    threshold = request.GET["threshold"]
    # filter_pa_pi = request.GET.get("checkbox") == "on" #TODO
    filter_pa_pi = True

    try:
        threshold = int(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")

    # query database
    targets = [urlparse.unquote(target) 
        for target in targets]

    target_hits, columns = get_target_hits(
        targets, threshold, filter_pa_pi=filter_pa_pi,
        as_dict=True)
    num_hits = len(target_hits)
    show_images = num_hits<MAX_HITS_FOR_IMAGE
    cols_to_remove = {"smiles"}
    if not show_images:
        cols_to_remove.add("image")
    columns = filter_columns(columns, cols_to_remove)

    request.session["targets"] = targets
    request.session["threshold"] = threshold
    request.session["hits"] = target_hits
    # request.session["columns"] = columns

    context = {
        "targets": targets,
        "threshold": threshold,
        "target_hits": target_hits,#[:MAX_RECORDS],
        "num_hits": num_hits,
        "columns": columns,
        "allow_optimise": True
    }

    return render(request,
        "natural_products/targets/target_hits.html", 
        context)