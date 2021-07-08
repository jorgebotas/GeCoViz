from django.http import HttpResponse, HttpResponseNotFound, JsonResponse
from django.conf import settings
from json import dumps

from .src.get_context import launch_analysis as gmgfam_query
from .src.gmgcFam_context import get_context as gmgcFam_query
from .src.gmgcFam_context import get_newick as gmgcFam_tree
from .src.query_gmgc import query_fam as query_gmgc
from .src.query_fam import get_neighborhood


RESULTS_PATH = settings.BASE_DIR + '/gmgfam/tmp/'

def get_context(request, datatype, query, cutoff):
    analysis = get_neighborhood("", query.split(','))
    return JsonResponse(analysis, safe=False)

def get_tree(request, query):
    try:
        with open(RESULTS_PATH + query + "_tree.nwx") as handle:
            newick = str(handle.read())
        return HttpResponse(newick, content_type='text/plain')
    except:
        print("NO TREE for specified query: " + str(query))
    return HttpResponseNotFound()

def get_eggNOG_levels(request):
    with open(RESULTS_PATH + 'eggNOG_LEVELS.txt') as handle:
        eggs = str(handle.read())
    return HttpResponse(eggs, content_type='text/plain')
