import numpy as np

from django.shortcuts import render

from django.http import HttpResponse

from .backend import query_pass, get_compound_info, draw_molecule

import html
import urllib.parse as urlparse

# Create your views here.

with open("PASS_targets.txt", "r") as f:
    pass_targets = [line.rstrip()
        for line in f.readlines()]

def index(request):
    context = {}
    return render(request, 
        "natural_products/index.html", context)
        
def target(request):
    targets = [
        (target,
            urlparse.quote(target))
        for target in pass_targets]
    thresholds = ["{:.02f}".format(threshold)
        for threshold in np.arange(0, 1, 0.05)[::-1]]
    context = {
        "targets": targets,
        "thresholds": thresholds}
    return render(request,
        "natural_products/targets.html", context)

def results(request, ):

    target = request.GET["target"]
    threshold = request.GET["threshold"]
    try:
        threshold = float(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")

    # query database
    target = urlparse.unquote(target)
    records = query_pass(target, threshold)

    context = {
        "target": target,
        "threshold": threshold,
        "records": records
    }

    return render(request,
        "natural_products/results.html", context)

def compound_info(request, compound_id):

    compound_id = "CNP" + compound_id

    # query database
    info, activities = get_compound_info(compound_id)

    context = {}

    assert isinstance(info, dict), type(info)
    assert "clean_smiles" in info
    print (info["clean_smiles"])

    smiles = info["clean_smiles"]
    if smiles is not None:
        img_filename = draw_molecule(smiles)
        context.update({"img_filename": img_filename})

    # convert to list of dicts for easy presentation
    info = [
        (k, v)
            for k, v in info.items()
    ]

    activities = [
        (target, activities[target])
            for target in pass_targets
    ]
    # sort activities by activity
    activities = sorted(activities, 
        key=lambda x: x[1], reverse=True)

    context.update({
        "info": info,
        "activities": activities})


    return render(request,
        "natural_products/compound_info.html", context)