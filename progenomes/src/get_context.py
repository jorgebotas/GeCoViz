from django.conf import settings
from json import dump
import pickle
from pymongo import MongoClient

from .Progenomes.main import make_tree as get_egg_context

STATIC_PATH = settings.BASE_DIR + '/static/geco/'
RESULTS_PATH = settings.BASE_DIR + '/progenomes/src/Progenomes/tmp/'

def get_pickle(filepath):
    """generate kegg pathway dictionary containning kegg descriptions"""
    # read dictionary from pickle file
    with open(filepath, 'rb') as pickle_in:
        pdict = pickle.load(pickle_in)
    return pdict

def get_egg_description(egg):
    """ connect to mongo progenomes5 database
    and retrieve description of Egnog """
    coll_e5 = MongoClient('10.0.3.1', 27017, maxPoolSize=10).gmgc_unigenes.progenomes_v5
    e5_desc = coll_e5.find({"e":egg})[0]["d"]
    return e5_desc

def get_info(l, egg_levels):
    gene = l[1]
    if gene == "_NA_":
        gene = "NA"
    info = { "gene" : gene.split(".")[-1],
              "start" : int(l[2]),
              "end" : int(l[3]),
              "size" : abs(int(l[3]) - int(l[2])),
              "strand" : l[4],
             }
    keggs = l[5].split("/")
    kegg = {}
    for k in keggs:
        if k == "_NA_":
            kegg = {}
            break
        description = ""
        try:
            description = KEGG_DICT[k]
        except: pass
        kegg[k] = { "description" : description }

    eggs = l[6].split("/")
    egg = {}
    for e in eggs:
        if e == "_NA_":
            egg = {}
            break
        i, l = e.split("@")
        try:
            lev = TAX_DICT[l]
            egg_levels[l] = lev
        except: pass
        try:
            description = get_egg_description(i)
        except:
            description = ""
        egg[l] = { i : { "id" : e,
                         "description" : description } }
    info["KEGG"] = kegg
    info["eggNOG"] = egg
    return info

def neigh2dict(csv_file, n_neigh):
    lines = []
    with open(csv_file, "r") as handler:
        for i, l in enumerate(handler.readlines()):
            if i != 0:
                lines.append(l.strip("\n").strip().split(","))
    data = {}
    neighborhood = {}
    central_gene = "NA"
    egg_levels = {}
    for i, l in enumerate(lines):
        pos = (i % (2 * n_neigh + 1)) - n_neigh
        neighborhood[pos] = get_info(l, egg_levels)
        if pos == n_neigh:
            central_gene = neighborhood[0]['gene']
            data[central_gene] = { "neighbourhood" : dict(neighborhood) }
            neighborhood = {}
            central_gene = "NA"
    return data, egg_levels

def get_context(query, n_neigh):
    get_egg_context(query)
    global KEGG_DICT, TAX_DICT
    kegg_path = STATIC_PATH + "pickle/KEGG_DESCRIPTION.pickle"
    KEGG_DICT = get_pickle(kegg_path)
    tax_path = STATIC_PATH + "pickle/TAX_LEVELS.pickle"
    TAX_DICT = get_pickle(tax_path)
    data, \
        egg_levels = neigh2dict(RESULTS_PATH + "genes_coordinates_sorted.csv", n_neigh)
    with open(RESULTS_PATH + 'eggNOG_LEVELS.txt', 'w') as handle:
        dump(egg_levels, handle)
    return dict(data)
