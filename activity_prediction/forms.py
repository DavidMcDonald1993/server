from django import forms 

class UploadFileForm(forms.Form):
    file_field = forms.FileField(
        label="Compounds (SDF or SMILES format)",
        max_length=1000)