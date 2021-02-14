from django.conf import settings
from django.shortcuts import render, redirect

from .forms import EggnogForm


def set_nightmode(request):
    request.session['nightmode'] = request.session.get('nightmode', False)
    mode = request.session['nightmode']
    if request.method == 'POST':
        request.POST['nightmode']
        request.session['nightmode'] = not mode
    return

def input(request):
    try:
        set_nightmode(request)
    except:
        if request.method == "POST":
            form = EggnogForm(request.POST)
            if form.is_valid():
                eggnog = form.cleaned_data['eggnog']
                return redirect('context',
                                query=eggnog)
    context = {
        'eggnog_form' : EggnogForm(),
    }
    return render(request, 'progenomes/input.html', context)

def context(request, query):
    try:
        set_nightmode(request)
    except: pass
    return render(request, 'progenomes/context.html', { 'query' : query })
