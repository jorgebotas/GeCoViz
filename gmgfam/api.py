from django.http import HttpResponse, HttpResponseNotFound, JsonResponse
from django.conf import settings
from json import dumps

from .src.get_context import launch_analysis as gmgfam_query
from .src.gmgcFam_context import get_context as gmgcFam_query


RESULTS_PATH = settings.BASE_DIR + '/gmgfam/results.tmp/'

# def get_context(request, datatype, query, cutoff):
    # if datatype == "cluster":
        # isCluster = True
    # else:
        # isCluster = False
    # if datatype == "list":
        # isList = True
        # query = query.split(',')
    # else:
        # isList = False
    # analysis_json = gmgfam_query(query,
                                # 20,
                                # cutoff,
                                # isCluster,
                                # isList)
    # return HttpResponse(analysis_json, content_type='application/json')

def get_context(request, datatype, query, cutoff):
    analysis = gmgcFam_query(query)
    return JsonResponse(analysis)


def get_tree(request, cluster):
    try:
        with open(RESULTS_PATH +
                  cluster + "_newick.txt") as handle:
            tree = str(handle.read())

        return HttpResponse(tree, content_type='text/plain')
    except:
        print("NO TREE for specified cluster: " + str(cluster))
    return HttpResponseNotFound()

def get_eggNOG_levels(request):
    with open(RESULTS_PATH + 'eggNOG_LEVELS.txt') as handle:
        eggs = str(handle.read())
    return HttpResponse(eggs, content_type='text/plain')
