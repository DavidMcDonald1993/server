from django import forms 

class UploadFileForm(forms.Form):
    file_field = forms.FileField(label="SDF_file:",
        max_length=100)