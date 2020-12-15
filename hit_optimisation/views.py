import os

from django.shortcuts import render

from django.views.static import serve

import multiprocessing as mp

from .forms import UploadFileForm
from .backend import hit_optimisation

# import pypdb
# from pypdb.clients.search.search_client import perform_search
# from pypdb.clients.search.search_client import SearchService, ReturnType
# from pypdb.clients.search.operators import text_operators

# search_service = SearchService.TEXT
# search_operator = text_operators.ExactMatchOperator(value="Homo sapiens",
#     attribute="rcsb_entity_source_organism.taxonomy_lineage.name")
# return_type = ReturnType.POLYMER_ENTITY

# human_targets = perform_search(search_service, search_operator, return_type)

# human_targets = {human_target[:4] 
#     for human_target in human_targets}

# Create your views here.
from django.http import HttpResponse, HttpResponseRedirect

def index(request):
    context = {}
    return render(request, 
        "hit_optimisation/index.html", context)

def upload(request):

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
        if form.is_valid():

            user_name = request.POST["user_name"]
            user_email = request.POST["user_email"]
            target = request.POST["target"]
            uploaded_file = request.FILES['file_field'] # name of attribute
            chain = request.POST["chain"]
            # num_generations = request.POST["num_generations"]
            user_settings = {key: int(request.POST[key]) 
                for key in settings}
           
            if uploaded_file.name.endswith(".smi"):
                # do optimisation
                # archive_filename = hit_optimisation(user_name, target, uploaded_file, chain, user_settings)
                # start new process that ends with sent email
                p  = mp.Process(target=hit_optimisation, args=(user_name, user_email, target, uploaded_file, chain, user_settings))
                p.start()
                print ("process spawned")
            #     return serve(request, 
            #         os.path.basename(archive_filename), 
            #         os.path.dirname(archive_filename))
                return HttpResponseRedirect("/hit_optimisation/success")
            else:
                form = UploadFileForm() # invalid sdf file
    else:
        form = UploadFileForm()

    context = {
        # "targets": human_targets,
        "form": form,
        "settings": settings.items()
    }

    if "target" in request.session:
        context.update({"target": request.session["target"]})

    if "smiles_filename" in request.session:
        context.update({"smiles_filename": request.session["smiles_filename"]})
    
    return render(request, 
        'hit_optimisation/upload.html', 
        context )

def success(request):

    context = {}

    return render(request,
        "hit_optimisation/success.html",
        context)