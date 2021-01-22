from django import forms 

class UploadFileForm(forms.Form):

    file_field = forms.FileField(label="Compound file (SDF or SMILES):", max_length=25)