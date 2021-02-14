from django.http import JsonResponse, HttpResponse, HttpResponseNotFound
from django.conf import settings

from ete3 import Tree
from json import dumps
import pickle

from .src.getContext import launch_analysis as eggnog_query


STATIC_PATH = settings.BASE_DIR + '/static/gecoviz/'
RESULTS_PATH = settings.BASE_DIR + '/progenomes/tmp/'

def get_context(request, query):
    analysis = eggnog_query(query, 2)
    # print(analysis)
    return JsonResponse(analysis, safe=False)

def get_pickle(filepath):
    """ read dictionary from pickle file """
    with open(filepath, 'rb') as pickle_in:
        pdict = pickle.load(pickle_in)
    return pdict

def get_tree(request, query):
    try:
        with open(RESULTS_PATH + query + "_tree.nwx") as handle:
            newick = str(handle.read())
        return HttpResponse(newick, content_type='text/plain')
    except:
        print("NO TREE for specified query: " + str(query))
    return HttpResponseNotFound()

def get_eggNOG_levels(request):
    return HttpResponseNotFound
    with open(RESULTS_PATH + 'eggNOG_LEVELS.txt') as handle:
        eggs = str(handle.read())
    return HttpResponse(eggs, content_type='text/plain')
