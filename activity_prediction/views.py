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

def index_view(request):
    context = {}
    return render(request, 
        "activity_prediction/index.html", context)

def login_page_view(request):
    if request.method == "POST":
        username = request.POST['username']
        password = request.POST['password']
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            # Redirect to a success page.
            return HttpResponseRedirect("/activity_prediction")
        else:
            # Return an 'invalid login' error message.
            return HttpResponseRedirect("/login_unsuccessful")
    else: 
        context = {}
        return render(request, "activity_prediction/login.html", context)

def logout_page_view(request):
    logout(request)
    return HttpResponseRedirect("/")

def login_unsuccessful_view(request):
    context = {}
    return render(request, 
    "activity_prediction/login_unsuccessful.html",
        context)

def upload_file_view(request):
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

            p = mp.Process(target=activity_predict,
                args=(request.user, uploaded_file),
                kwargs={"threshold": threshold})
            p.start()
            print ("process spawned")

            return HttpResponseRedirect("/activity_prediction/success")
                
    else:
        form = UploadFileForm()

    context = {"form": form}
    context["username"] = request.user.username
    context["user_email"] = request.user.email
    
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
            # return HttpResponseRedirect("/")
            # if "authenticated" in request.session.keys() and request.session["authenticated"]:
            #     del request.session["authenticated"]
            #     filename = get_file_from_token(token, user.id)
            #     if filename is None:
            #         return HttpResponseRedirect("/download_error")
            #     response = FileResponse(open(filename, 'rb'))
            #     return response

            context["authenticated"] = True # change page 
            request.session["authenticated_user_id"] = user.id

            
            # context["filename"] = filename
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
