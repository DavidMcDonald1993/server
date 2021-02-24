import os

from django.shortcuts import render

from django.core.exceptions import ObjectDoesNotExist

# Create your views here.
from django.http import HttpResponse, HttpResponseRedirect, FileResponse
from django.views.static import serve

from activity_prediction.forms import UploadFileForm
from activity_prediction.backend import activity_predict

from utils.users import get_file_from_token

import multiprocessing as mp

from django.contrib.auth import authenticate, login, logout

import pandas as pd
import zipfile

# from utils.io import load_json
import json

def index_view(request):
    if not request.user.is_authenticated:
        return HttpResponseRedirect("/login")
    context = {}
    return render(request, 
        "activity_prediction/index.html", context)

def login_view(request):
    if request.method == "POST":
        username = request.POST['username']
        password = request.POST['password']
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            # Redirect to a success page.
            return HttpResponseRedirect("/")
        else:
            # Return an 'invalid login' error message.
            return HttpResponseRedirect("/login_unsuccessful")
    else: 
        context = {}
        return render(request, "activity_prediction/login.html", context)

def logout_view(request):
    logout(request)
    return HttpResponseRedirect("/")

def login_unsuccessful_view(request):
    context = {}
    return render(request, 
    "activity_prediction/login_unsuccessful.html",
        context)

def upload_file_view(request):
    if not request.user.is_authenticated:
        return HttpResponseRedirect("/login")

    ppb2_options = [
        ("Algorithm 1", "morg2-nn+nb"),
        ("Algorithm 2", "morg3-xgc"),
    ]

    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            uploaded_file = request.FILES["file_field"] # name of attribute
            threshold = int(request.POST["threshold"])
            use_pass = request.POST.get("use_pass") == "on"
            # use_pass = True
            use_ppb = request.POST.get("use_ppb") == "on"

            if use_ppb:
                model = request.POST["ppb2_option"]
            else:
                model = None

            # use_ppb = False
            perform_enrichment = (use_pass or use_ppb) and request.POST.get("perform_enrichment") == "on"
            group_compounds = request.POST.get("group_compounds") =="on"

            # handle with multi processing 
            p = mp.Process(target=activity_predict,
                args=(request.user, uploaded_file),
                kwargs={
                    "enrichment_threshold": threshold, 
                    "pass_predict": use_pass, 
                    "ppb2_predict": use_ppb,
                    "model": model,
                    "perform_enrichment": perform_enrichment,
                    "group_compounds": group_compounds
                })
            p.start()
            print ("process spawned")

            return HttpResponseRedirect("/activity_prediction/success")
                
    else:
        form = UploadFileForm()

    context = {
        "form": form,
        "username": request.user.username,
        "user_email": request.user.email,
        "model_choices": ppb2_options
    }
    
    return render(request, 
        'activity_prediction/upload.html', 
        context)

def success_view(request):

    context = {}

    return render(request, 
        "activity_prediction/success.html",
        context)

def download_view(request, token):

    context = {"token": token}

    if request.method == "POST":

        username = request.POST['username']
        password = request.POST['password']
        user = authenticate(request, username=username, password=password)

        if user is not None:

            context["authenticated"] = True # change page 
            request.session["authenticated_user_id"] = user.id

            # display summary on download screen
            filename = get_file_from_token(token, user.id)

            if "activity_prediction" in filename:
                zf = zipfile.ZipFile(filename)

                # all compounds in 

                compounds = sorted(filter(lambda s: s!="", 
                    set([os.path.dirname(f) for f in zf.namelist() ])))
                if "compound_group" in compounds:
                    compounds = ["compound_group"]

                summary = {
                    compound: dict() for compound in compounds
                }

                # predicted targets
                # context["predicted_targets"] = json.loads(
                #     zf.open("uniprot_confidences.json").read())
                for compound in compounds:
                    summary[compound]["Targets"] = \
                        pd.read_csv(zf.open(f"{compound}/combined_uniprot_confidences.tsv"),
                            sep="\t")
                    summary[compound]["Enrichment"] = \
                        pd.read_csv(zf.open(f"{compound}/enrichment.csv"),
                            sep=",")
                    summary[compound]["SimilarDrugs"] = \
                        pd.read_csv(zf.open(f"{compound}/similar_drugs.tsv"),
                            sep="\t")
                    summary[compound]["AssociatedDisease"] = \
                        pd.read_csv(zf.open(f"{compound}/associated_diseases.tsv"),
                            sep="\t")
                context["summary"] = summary

                # # enrichment dfs
                # context["enrichment_dfs"] = []
                # for filename in filter(lambda z: 
                #         z.filename.endswith("enrichment.csv"), zf.filelist):
                #     df = pd.read_csv(zf.open(filename.filename))
                #     context["enrichment_dfs"].append(
                #         (filename.filename, df.head()))

        else:
            context["login_error"] = True
    
    elif "authenticated_user_id" in request.session.keys():
        authenticated_user_id = request.session["authenticated_user_id"]
        del request.session["authenticated_user_id"]
        filename = get_file_from_token(token, authenticated_user_id)
        if filename is None:
            return HttpResponseRedirect("/download_error")
        response = FileResponse(open(filename, 'rb'))
        return response

    return render(request, 
        "activity_prediction/download.html",
        context)

def download_error_view(request):
    
    context = {}

    return render(request, 
        "activity_prediction/download_error.html",
        context)
