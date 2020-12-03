import numpy as np

from django.shortcuts import render

from django.http import HttpResponse

from .backend import query_pass_activities, get_compound_info, draw_molecule, get_multiple_compound_info

import html
import urllib.parse as urlparse

from .backend import get_categories, get_targets_for_category

# Create your views here.

# with open("PASS_targets.txt", "r") as f:
#     pass_targets = [line.rstrip()
#         for line in f.readlines()]

def index(request):
    context = {}
    return render(request, 
        "natural_products/index.html", context)
        
def target(request):
    # targets = [
    #     (target,
    #         urlparse.quote(target))
    #     for target in pass_targets]

    categories = get_categories()
    targets = [(category, 
            [(target, urlparse.quote(target))
                for target in get_targets_for_category(category)])
        for category in categories]

    thresholds = ["{:.02f}".format(threshold)
        for threshold in np.arange(0, 1, 0.05)[::-1]]
    context = {
        "targets": targets,
        "thresholds": thresholds}
    return render(request,
        "natural_products/targets.html", context)

def results(request, ):

    category = request.GET["category"]
    target = request.GET[category] # get values of fields by id
    threshold = request.GET["threshold"]
    filter_pa_pi = request.GET.get("checkbox") == "on"

    try:
        threshold = float(threshold)
    except ValueError:
        return HttpResponse("Invalid threshold")

    # query database
    target = urlparse.unquote(target)
    records = query_pass_activities(category, target, threshold, filter_pa_pi=filter_pa_pi)

    context = {
        "target": target,
        "threshold": threshold,
        "records": records,
        "num_hits": len(records)
    }

    return render(request,
        "natural_products/results.html", context)

def all_compounds(request):

    compounds = get_multiple_compound_info()

    # get compounds molecular formulas from database
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

    smiles = info["clean_smiles"]
    if smiles is not None:
        img_filename = draw_molecule(smiles)
        if img_filename is not None:
            context.update({"img_filename": img_filename})

    # convert to list of dicts for easy presentation
    info = [
        (k, v)
            for k, v in info.items()
    ]

    # if activities is not None:
    #     activities = [
    #         (target, activities[target])
    #             for target in pass_targets
    #     ]
    #     # sort activities by activity
    #     activities = sorted(activities, 
    #         key=lambda x: x[1], reverse=True)

    context.update({
        "compound_id": compound_id,
        "info": info,
        "activities": activities})

    return render(request,
        "natural_products/compound_info.html", context)