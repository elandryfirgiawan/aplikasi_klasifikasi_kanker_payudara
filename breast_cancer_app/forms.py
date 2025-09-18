from django import forms

class UploadFileForm(forms.Form):
    file = forms.FileField(label="Upload File Excel")
    
class UploadFastaForm(forms.Form):
    fasta_file = forms.FileField(label="Upload File FASTA")
