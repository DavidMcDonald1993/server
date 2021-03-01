import os

import numpy as np
import pandas as pd

from django.shortcuts import render

from django.http import HttpResponse, HttpResponseRedirect
from django.views.static import serve

import html
import urllib.parse as urlparse

from natural_products.backend import (
    get_coconut_compound_info_from_mongo,
)
from utils.queries import (
    get_all_activities_for_compound,
    get_info_for_multiple_compounds,
    get_combined_uniprot_confidences_for_compounds,
    get_drugs_for_uniprots,
    get_diseases_for_uniprots,
    get_all_pathways_for_uniprots,
    get_all_reactions_for_uniprots,
)

from natural_products.views import DEFAULT_THRESHOLD, MAX_HITS_FOR_IMAGE

# Create your views here.
def all_compounds_view(request):

    context = {}

    if request.method == 'POST':

        name_like = request.POST["name_like"]
        formula_like = request.POST["formula_like"]
        smiles_like = request.POST["smiles_like"]

        if name_like == "":
            name_like = None 
        if formula_like == "":
            formula_like = None
        if smiles_like == "":
            smiles_like = None

        columns = [
            "image_path",
            "coconut_id",
        ]

        show_name = request.POST.get("show_name") == "on"
        show_formula = request.POST.get("show_formula") == "on"
        show_smiles = request.POST.get("show_smiles") == "on"

        if show_name:
            columns.append("name")
        if show_formula:
            columns.append("formula")
        if show_smiles:
            columns.append("smiles")

        compounds = get_info_for_multiple_compounds(
            name_like=name_like,
            formula_like=formula_like,
            smiles_like=smiles_like,
            columns=columns,
        )

        num_hits = len(compounds)

        compounds = [
            {k: v for k, v in zip(columns, hit)}
            for hit in compounds
        ]

        context["show_results"] = True
        context["compounds"] = compounds
        context["show_name"] = show_name
        context["show_formula"] = show_formula
        context["show_smiles"] = show_smiles
        context["show_images"] = num_hits < MAX_HITS_FOR_IMAGE

    return render(request, 
        "natural_products/compounds/all_compounds.html", context)

def compound_info_view(request, compound_id):

    compound_id = "CNP" + compound_id
    threshold = DEFAULT_THRESHOLD

    context = {
        "compound_id": compound_id,
        "threshold": threshold,
    }

    # query database
    compound_info = get_coconut_compound_info_from_mongo(compound_id)

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

    uniprots, uniprot_cols = get_combined_uniprot_confidences_for_compounds(compound_id,
        threshold=threshold, as_dict=True)

    all_accs = {
        record["acc"]: (record["confidence_score"], record["confidence_type"])
        for record in uniprots
    }

    if len(all_accs) > 0:
        all_acc_list = list(all_accs)
        pathways = get_all_pathways_for_uniprots(all_acc_list)
        reactions = get_all_reactions_for_uniprots(all_acc_list)
        '''
        add confidences to drugs and diseases
        '''
        drugs, drug_cols = get_drugs_for_uniprots(all_acc_list, as_dict=False)
        '''
        d.drug_name, d.inchi, 
        d.canonical_smiles, d.drug_type, d.drug_class, d.company, 
        disease.disease_name, drd.clinical_status,
        u.acc, 
        ud.activity, 
        ud.reference
        '''
        drugs = [
            (drug_name, inchi, smiles, drug_type, drug_class, company,
                disease_name, status, acc, f"{all_accs[acc][0]}", f"{all_accs[acc][1]}", activity, reference)
            for drug_name, inchi, smiles, drug_type, drug_class, company,\
                disease_name, status, acc, activity, reference in drugs
        ]
        diseases, disease_cols = get_diseases_for_uniprots(all_acc_list, as_dict=False)
        '''
        d.disease_name, d.icd,
        u.acc, ud.clinical_status
        '''
        diseases = [
            (disease_name, icd, acc, f"{all_accs[acc][0]}", f"{all_accs[acc][1]}", status)
            for disease_name, icd, acc, status in diseases
        ]
    else:
        pathways = []
        reactions = []
        drugs = []
        diseases = []

    context.update({
        "compound_name": compound_name,
        "compound_info": compound_info,
        "activities": activities,
        # "inferred_uniprots": inferred_uniprots,
        # "predicted_uniprots": predicted_uniprots,
        "uniprots": uniprots,
        "pathways": pathways,
        "reactions": reactions,
        "drugs": drugs,
        "diseases": diseases,
    })

    return render(request,
        "natural_products/compounds/compound_info.html", 
        context)






