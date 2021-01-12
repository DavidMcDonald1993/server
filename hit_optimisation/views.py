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


# human_targets = get_human_targets()

# Create your views here.
from django.http import HttpResponse, HttpResponseRedirect

def index_view(request):
    context = {}
    return render(request, 
        "hit_optimisation/index.html", context)

def upload_view(request):

    assert request.user.is_authenticated

    settings = {
        "number_of_mutants_first_generation": 10,
        "number_of_crossovers_first_generation": 10,
        "number_of_mutants": 3,
        "number_of_crossovers": 3,
        "number_elitism_advance_from_previous_gen": 3,
        "top_mols_to_seed_next_generation": 5,
        "diversity_mols_to_seed_first_generation": 3,
        "diversity_seed_depreciation_per_gen": 0,
        "num_generations": 100,
    }

    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if "smiles_filename" in request.session or form.is_valid():

            target = request.POST["target"]
            if "smiles_filename" in request.session:
                uploaded_file = request.session["smiles_filename"] # string full location on server
                del request.session["smiles_filename"]
            else:
                # smiles file from client
                uploaded_file = request.FILES["file_field"] # name of attribute
           
            chain = request.POST["chain"]
           
            user_settings = {key: int(request.POST[key]) for key in settings}
           
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
        "settings": settings.items(),
        "form": form,
        "username": request.user.username,
        "user_email": request.user.email,
    }

    if "targets" in request.session:
        targets = request.session["targets"]
        targets_to_gene_symbols = search_for_targets(targets, only_best_scoring=False)

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