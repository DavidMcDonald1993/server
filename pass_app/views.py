import os

from django.shortcuts import render

from django.core.exceptions import ObjectDoesNotExist

# Create your views here.
from django.http import HttpResponse, HttpResponseRedirect, FileResponse
from django.views.static import serve

from pass_app.forms import UploadFileForm
from pass_app.backend import pass_predict

from utils.users import get_file_from_token

import multiprocessing as mp

from django.contrib.auth import authenticate, login, logout

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

def logout_page(request):
    logout(request)
    return HttpResponseRedirect("/")

def login_unsuccessful(request):
    context = {}
    return render(request, 
    "pass_app/login_unsuccessful.html",
        context)

def upload_file(request):
    if not request.user.is_authenticated:
        return HttpResponseRedirect("login_unsucessful")

    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            # username = request.POST["username"]
            # user_email = request.POST["user_email"]
            uploaded_file = request.FILES["file_field"] # name of attribute
            threshold = int(request.POST["threshold"])

            # handle with multi processing 

            p = mp.Process(target=pass_predict,
                args=(request.user, uploaded_file),
                kwargs={"threshold": threshold})
            p.start()
            print ("process spawned")

            return HttpResponseRedirect("/pass_app/success")
                
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

def download(request, token):

    context = {"token": token}

    if request.method == "POST":

        username = request.POST['username']
        password = request.POST['password']
        user = authenticate(request, username=username, password=password)

        if user is not None:
            # return HttpResponseRedirect("/")
            filename = get_file_from_token(token, user.id)
            if filename is None:
                return HttpResponseRedirect("/download_error")
            # if filename is not None:
            # response = FileResponse(open(filename, 'rb'))
            # return response
            context["filename"] = filename
        else:
            context["login_error"] = True

    return render(request, 
        "pass_app/download.html",
        context)

def download_error(request):
    
    context = {}

    return render(request, 
        "pass_app/download_error.html",
        context)

def favicon(request):
    return HttpResponse("/favicon")
