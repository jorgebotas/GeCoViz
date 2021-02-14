from django.conf import settings
from django.shortcuts import render, redirect

from .forms import ClusterForm, ListForm, NovelFamForm

def set_nightmode(request):
    request.session['nightmode'] = request.session.get('nightmode', True)
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
            if request.POST['input_choice'] == "cluster":
                # form = ClusterForm(request.POST)
                # if form.is_valid():
                    # cluster = form.cleaned_data['unigene_cluster']
                    # cutoff = form.cleaned_data['cutoff']
                    # if cluster != "" and cutoff != "":
                        # return redirect('cluster_context',
                                        # cluster=cluster,
                                        # cutoff=cutoff)
                form = NovelFamForm(request.POST)
                if form.is_valid():
                    cluster = form.cleaned_data['unigene_cluster']
                    cutoff = 30
                    if cluster != "" and cutoff != "":
                        return redirect('cluster_context',
                                        cluster=cluster,
                                        cutoff=cutoff)

            elif request.POST['input_choice'] == "list":
                form = ListForm(request.POST, request.FILES)
                if form.is_valid():
                    gene_list = form.cleaned_data['gene_list'].strip().split('\n')
                    if "".join(gene_list).strip() != "":
                        gene_list = [g.strip("\n").strip("\t").strip("") for g in
                                     gene_list]
                        genelist = ",".join(gene_list)
                        return redirect('list_context', genelist=genelist)
    context = {
        # 'cluster_form' : ClusterForm(),
        # 'list_form' : ListForm(),
        'novelFam_form' : NovelFamForm(),
    }
    return render(request, 'gmgfam/input.html', context)

def cluster_context(request, cluster, cutoff):
    try:
        set_nightmode(request)
    except: pass
    return render(request, 'gmgfam/context.html', { 'query' : cluster,
                                                  'cutoff' : cutoff,
                                                  'data' : ""
                                                 })

def unigene_context(request, unigene):
    try:
        set_nightmode(request)
    except:
        pass
    return render(request, 'gmgfam/context.html', { 'query' : unigene,
                                                  'cutoff' : 0,
                                                  'data' : ""
                                                 })

def list_context(request, genelist):
    try:
        set_nightmode(request)
    except:
        pass
    return render(request, 'gmgfam/context.html', {
                                            'genelist' : genelist,
                                                 })
