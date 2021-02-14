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
    return JsonResponse(analysis)

def get_pickle(filepath):
    """ read dictionary from pickle file """
    with open(filepath, 'rb') as pickle_in:
        pdict = pickle.load(pickle_in)
    return pdict

def get_tree(request, query):
    try:
        with open(RESULTS_PATH + "tree.nwx") as handle:
            tree = str(handle.read())
        t = Tree(tree)
        tax_path = STATIC_PATH + "pickle/TAX_LEVELS.pickle"
        tax_dict = get_pickle(tax_path)
        print(tax_dict['1224318'])
        for node in t:
            name = str(node.name).split(".")
            try:
                tax = tax_dict[name[0]]
                tax = tax.replace(".", "").strip().split(" ")
                tax_s = "_".join(tax[:2])
                tax_l = "_".join(tax)
                node.name = ".".join([tax_s] + name[:1] + [tax_l] + name[1:])
            except: pass
        newick = t.write()
        return HttpResponse(newick, content_type='text/plain')
    except:
        print("NO TREE for specified query: " + str(query))
    return HttpResponseNotFound()

def get_eggNOG_levels(request):
    with open(RESULTS_PATH + 'eggNOG_LEVELS.txt') as handle:
        eggs = str(handle.read())
    return HttpResponse(eggs, content_type='text/plain')
