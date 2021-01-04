from django import forms 

class UploadFileForm(forms.Form):

    # def __init__(self):
        # super(UploadFileForm, self).__init__()
    file_field = forms.FileField(label="smiles_file:",
        max_length=100)
        # self.file_field = file_field