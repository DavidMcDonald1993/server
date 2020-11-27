from django.shortcuts import render

from .forms import UploadFileForm
from .backend import hit_optimisation

# Create your views here.

def index(request):
    context = {}
    return render(request, 
        "hit_optimisation/index.html", context)

def upload(request):

    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():

            target = request.POST["target"]
            uploaded_file = request.FILES['file_field'] # name of attribute
           
            if uploaded_file.name.endswith(".smi"):
                # do optimisation
                optimised_hits = hit_optimisation(target, uploaded_file)
            else:
                form = UploadFileForm() # invalid sdf file
    else:
        form = UploadFileForm()
    return render(request, 
        'hit_optimisation/upload.html', 
        {'form': form})
