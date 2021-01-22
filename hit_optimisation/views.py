import os

from django.shortcuts import render

from django.views.static import serve

import multiprocessing as mp

import urllib.parse as urlparse


from .forms import UploadFileForm
from .backend import hit_optimisation

from utils.pdb_utils import get_pdb_ids_from_gene_symbol, get_human_targets
from utils.genenames_utils import search_for_targets

import urllib.parse as urlparse

from hit_optimisation.backend import load_base_settings, load_parameter_descriptions

# human_targets = get_human_targets()

# Create your views here.
from django.http import HttpResponse, HttpResponseRedirect

display_keys = [
    "number_of_mutants_first_generation",
    "number_of_crossovers_first_generation",
    "number_of_mutants",
    "number_of_crossovers",
    "number_elitism_advance_from_previous_gen",
    "top_mols_to_seed_next_generation",
    "diversity_mols_to_seed_first_generation",
    "diversity_seed_depreciation_per_gen",
    "num_generations",
    "max_variants_per_compound",
]

def upload_view(request):

    if not request.user.is_authenticated:
        return HttpResponseRedirect("/login")

    base_settings = load_base_settings()
    parameter_descriptions = load_parameter_descriptions()

    display_settings = [
        (key, parameter_descriptions[key], base_settings[key])
        for key in display_keys
    ]

    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if "smiles_filename" in request.session.keys() or form.is_valid():

            target = request.POST["target"]
            if "smiles_filename" in request.session.keys():
                uploaded_file = request.session.pop("smiles_filename") # string full location on server
            else:
                # smiles file from client
                uploaded_file = request.FILES["file_field"] # name of attribute
           
            chain = request.POST["chain"]
           
            user_settings = {key: int(request.POST[key]) 
                for key, _, _ in display_settings}
           
            if isinstance(uploaded_file, str) and uploaded_file.endswith(".smi")\
                    or uploaded_file.name.endswith(".smi"):
                # do optimisation
                # start new process that ends with sent email
                p  = mp.Process(target=hit_optimisation, 
                    args=(request.user, target, uploaded_file, chain, user_settings))
                p.start()
                print ("process spawned")
                return HttpResponseRedirect("/hit_optimisation/success")
            else:
                form = UploadFileForm() # invalid sdf file
    else:
        form = UploadFileForm()

    context = {
        # "targets": human_targets,
        "settings": display_settings,
        "form": form,
        "username": request.user.username,
        "user_email": request.user.email,
    }

    if "optimise" in request.session.keys() and request.session.pop("optimise") and "targets" in request.session.keys():
        targets = request.session.pop("targets")
        targets_to_gene_symbols = search_for_targets(targets, 
            only_best_scoring=True, max_hits=100)

        pdb_ids = sorted(
            ((target, symbol, score, pdb_id)
                for target in targets_to_gene_symbols
                for symbol, score in targets_to_gene_symbols[target]
                for pdb_id in get_pdb_ids_from_gene_symbol(symbol)),
        key=lambda t: t[3]) # sort by pdb_ID

        # pdb_ids = {pdb_id for gene_symbol in gene_symbols
            # for pdb_id in get_pdb_ids_from_gene_symbol(gene_symbol)}
        if len(pdb_ids) > 0:
            context["pdb_ids"] = pdb_ids#.intersection(human_targets)

    if "smiles_filename" in request.session:
        context["smiles_filename"] = \
            os.path.basename(request.session["smiles_filename"])

    return render(request, 
        'hit_optimisation/upload.html', 
        context )

def success_view(request):

    context = {}

    return render(request,
        "hit_optimisation/success.html",
        context)