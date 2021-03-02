from django import forms
from django.core.validators import MinValueValidator, MaxValueValidator

class ClusterForm(forms.Form):
    unigene_cluster = forms.CharField(max_length=11,
                                      label="",
                                      label_suffix="",
                                      required=False,
                                      widget=forms.TextInput(attrs={
                                          'placeholder' : 'XXX_XXX_XXX',
                                          'class' : 'form-control'
                                      }))
    cutoff = forms.FloatField(label="Cutoff percentage (default: 30%)",
                              label_suffix="",
                              required=False,
                              validators=[MinValueValidator(0),
                                          MaxValueValidator(100)],
                              widget=forms.NumberInput(attrs={
                                  'pattern' : '[0-100]',
                                  'value' : '30',
                                  'class' : 'form-control'
                              }))

class ListForm(forms.Form):
    gene_list = forms.CharField(label="",
                                label_suffix="",
                                required=False,
                                widget=forms.Textarea(attrs={
                                     'class' : 'form-control'
                                 }))

class FileForm(forms.Form):
    input_file = forms.FileField(label="",
                                 label_suffix="",
                                 required=False,
                                 widget=forms.FileInput(attrs={
                                     'accept' : '.json',
                                     'id' : 'input_file'
                                 }))

class NewickForm(forms.Form):
    newick_file = forms.FileField(label="",
                                 label_suffix="",
                                 required=False,
                                 widget=forms.FileInput(attrs={
                                     'accept' : '.nwx',
                                     'id' : 'input_newick'
                                 }))
