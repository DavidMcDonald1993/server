import os

from django.shortcuts import render

from django.core.exceptions import ObjectDoesNotExist

# Create your views here.
from django.http import HttpResponse, HttpResponseRedirect
from django.views.static import serve


from .models import Compound
from .forms import UploadFileForm
from .backend import pass_predict

import multiprocessing as mp

def index(request):
    # latest_compound_list = Compound.objects\
    #     .order_by("-added_date")[:5]
    # context = {"latest_compound_list": latest_compound_list}
    context = {}
    return render(request, 
        "pass_app/index.html", context)

def compound_detail(request, compound_name):
    try:
        compound = Compound.objects.get(name=compound_name)
        return HttpResponse(
            "you are looking at details " +\
            "for compound: " + compound_name +\
            "<br> smiles: " + compound.canonical_smiles)
    except ObjectDoesNotExist:
        return HttpResponse(compound_name +\
            " does not exist in database")

def upload_file(request):
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            user_name = request.POST["user_name"]
            user_email = request.POST["user_email"]
            uploaded_file = request.FILES['file_field'] # name of attribute
            if uploaded_file.name.endswith(".sdf"):
                # filepath = handle_uploaded_file(uploaded_file)
                # return serve(request, 
                #     os.path.basename(filepath), 
                #     os.path.dirname(filepath))

                # handle with multi processing 

                p = mp.Process(target=pass_predict,
                    args=(user_name, user_email, uploaded_file))
                p.start()

                print ("process spawned")

                return HttpResponseRedirect("/pass_app/success")
                
            else:
                form = UploadFileForm() # invalid sdf file
    else:
        form = UploadFileForm()
    return render(request, 
        'pass_app/upload.html', 
        {'form': form})

def success(request):

    context = {}

    return render(request, 
        "pass_app/success.html",
        context)

def favicon(request):
    return HttpResponse("favicon")
