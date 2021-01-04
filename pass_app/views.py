import os

from django.shortcuts import render

from django.core.exceptions import ObjectDoesNotExist

# Create your views here.
from django.http import HttpResponse, HttpResponseRedirect
from django.views.static import serve

from pass_app.forms import UploadFileForm
from pass_app.backend import pass_predict

import multiprocessing as mp

from django.contrib.auth import authenticate, login

def index(request):
    context = {}
    return render(request, 
        "pass_app/index.html", context)

def login_page(request):
    if request.method == "POST":
        username = request.POST['username']
        password = request.POST['password']
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            # Redirect to a success page.
            return HttpResponseRedirect("/pass_app")
        else:
            # Return an 'invalid login' error message.
            return HttpResponseRedirect("/login_unsuccessful")
    else: 
        context = {}
        return render(request, "pass_app/login.html", context)

def login_unsuccessful(request):
    context = {}
    return render(request, 
    "pass_app/login_unsuccessful.html",
        context)

def upload_file(request):
    assert request.user.is_authenticated
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            user_name = request.POST["user_name"]
            user_email = request.POST["user_email"]
            if "user_name" not in request.session:
                request.session["user_name"] = user_name
            if "user_email" not in request.session:
                request.session["user_email"] = user_email
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

    context = {"form": form}

    context["username"] = request.user.username
    context["user_email"] = request.user.email
    
    return render(request, 
        'pass_app/upload.html', 
        context)

def success(request):

    context = {}

    return render(request, 
        "pass_app/success.html",
        context)

def favicon(request):
    return HttpResponse("favicon")
