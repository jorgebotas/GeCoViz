from django import forms
from django.core.validators import MinValueValidator, MaxValueValidator

class EggnogForm(forms.Form):
    eggnog = forms.CharField(max_length=7,
                              label="",
                              label_suffix="",
                              required=True,
                              widget=forms.TextInput(attrs={
                                  'placeholder' : 'COGXXXX | XXXXX',
                                  'class' : 'form-control',
                              }))
