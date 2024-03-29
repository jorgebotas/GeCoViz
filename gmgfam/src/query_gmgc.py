from collections import defaultdict
from ete3 import Tree, NCBITaxa
import gridfs
from os import makedirs, path
import pickle
from pymongo import MongoClient
from json import load, dump
from shutil import rmtree
import time

from django.conf import settings

RESULTS_PATH = settings.BASE_DIR + '/gmgfam/tmp/'
STATIC_PATH = settings.BASE_DIR + '/static/gecoviz/'

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', 
        decimals = 1, length = 100, fill ='#', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()


def mongo_connect():
    """
    Connection to MongoDB
    """
    db = None
    if not db:
        client = MongoClient('10.0.3.1', 27017, maxPoolSize=10)
        db = client.gmgc_unigenes
        coll_unigene = db.neighbour
        coll_clusters = db.emapper_v2
        coll_e5 = db.eggnog_v5
        coll_taxa = db.taxo_map
        coll_seqs = db.sequences

    return [client, db, coll_unigene, coll_clusters, 
            coll_e5, coll_taxa, coll_seqs]


def get_pickle(filepath):
    """generate kegg pathway dictionary containning kegg descriptions"""
    # read dictionary from pickle file
    with open(filepath, 'rb') as pickle_in:
        pdict = pickle.load(pickle_in)
    return pdict


def clean_unigene(gmgc):
        '''in case GMGC nomeclature was used in the input example,
         GMGC.100_000_123.UNKNOWN words after/before . are removed'''
        gmgc_clean = gmgc.split(".")
        if len(gmgc_clean) > 1:
            return gmgc_clean[1]
        return gmgc


def cl_to_unigene(cl):
    u = coll_clusters.find_one({'cl' : cl}) or {}
    return u.get('u', cl)


def unigene_to_cl(unigene):
    cl = coll_clusters.find_one({'u' : unigene}) or {}
    return cl.get('cl', unigene)


def get_taxonomic_prediction(unigenes):
    tax_pred = defaultdict(list)

    tax_ids = coll_taxa.find({'u': {'$in': unigenes}}, {'_id': 0})
    for t in tax_ids:
        unigene = t['u']
        txid = t['txid']
        lineage = ncbi.get_lineage(int(txid))
        for tid in lineage:
            rank = ncbi.get_rank([tid]).values()
            rank = str(rank)[14:-3]
            if rank == 'no rank': 
                continue
            desc = ncbi.get_taxid_translator([tid]).values()
            desc = str(desc)[14:-3]

            tax_pred[unigene].append({
                'id': tid,
                'level': rank,
                'description': desc,
            })

    return tax_pred


def get_og_description(ogs):
    desc = coll_e5.find({'e': {'$in': ogs}})
    return {d['e']: d['d'] for d in desc}


def get_functional_info(unigenes):
    clusters = coll_clusters.find({'u': {'$in': unigenes}}, {'_id': 0})
    clusters = { c['u']: c for c in clusters }

    all_ogs = set() # Save all OGs to later get their description
    info = {}

    for unigene, cl in clusters.items():
        kpaths = []
        for k in cl.get('K_P', '').split(','):
            k_desc = kegg_dict.get(k)
            if k_desc:
                kpaths.append({'id': k, 'description': k_desc})
        ogs = []
        for og_lev in cl.get('OGs', '').split(','):
            og, level = og_lev.split('@')
            level = int(level)
            if level != 1:
                all_ogs.add(og)
                ogs.append({
                    'id': og,
                    'level': level,
                })
        info[unigene] = {
                'Gene name': cl.get('p_n', ''),
                'Orthologous groups': ogs,
                'KEGG pathway': kpaths,
                'GMGFam': cl.get('cl', ''),
        }

    # Add eggNOG orthologous group description
    og_description = get_og_description(list(all_ogs))
    for i in info.values():
        for og in i['Orthologous groups']:
            og['description'] = og_description.get(og['id'], '')

    return info


def write_newick(query_cluster, results_dir):
    treedb = client.trees
    fs = gridfs.GridFS(treedb)
    tfile = fs.find_one({"filename":query_cluster})
    if not tfile:
        return
    treedata = fs.get(tfile._id)
    orig_newick = str(treedata.read(), "utf-8:")
    # Clean up newick
    t = Tree(orig_newick)
    cluster_unigenes = []
    for node in t:
        node.name = node.name.split('.')[1]
        cluster_unigenes.append(node.name)
    newick = t.write()
    with open(results_dir + query_cluster + "_tree.nwx", "w") as outputfile:
        outputfile.write(newick)
    return cluster_unigenes


def swap_strand(s, reference_s):
    if reference_s == "+":
        return s
    else:
        if s == "+":
            return "-"
        elif s == "-":
            return "+"
    return "+"


def swap_strands(gene_info, reference_s):
    for gene in gene_info.values():
        gene['strand'] = swap_strand(gene['strand'], reference_s)


def get_members(cl):
    members = client.gmgc_clusters.members.find_one({'cl' : cl}) or {}
    return members.get('clm', [])


def get_orfs(unigenes):
    """ Retrieve orf genomic information
    Returns dict(unigene: dict(central_gene: info))
    """
    cluster = coll_unigenes.find({ 'u': {'$in': unigenes} })
    orfs = defaultdict(dict)
    for orf in cluster:
        unigene = orf['u']
        orf = orf['o']
        for a in orf[:max_orfs]:
            gene = a['g']
            locus = a['s']
            start = locus[0]
            end = locus[1]
            strand = locus[2]
            try:
                orfs[unigene][gene]=[start,end,strand]
            except:
                print("ERROR retrieving ORF info from gmgc cluster")
    return orfs


def get_neighbors(orfs):
    """ get a dict containing -2-1+1+2 neighbor genes list surrounding
    the unigenes in every contig, corrected by strand orientation """
    neigh_dict = {}
    for k,v in orfs.items():
        gene_ordered =[]
        query_gene = k
        start = v[0]
        end = v[1]
        strand = v[2]
        orf = k.split("_")
        gene = int(orf[3])
        sample_contig = orf[0:3]

        for gene_pos in range(-neighbor_range,neighbor_range+1,1):
                sample_cont = []
                genes = int(gene)+int(gene_pos)
                sample_cont = sample_contig[0:3]
                sample_cont.append(str(genes))
                genes = "_".join(sample_cont)
                gene_ordered.append(genes)
        neigh_dict[query_gene] = {
            'strand' : strand,
            'neighborhood' : gene_ordered
        }
    return neigh_dict


def get_gene_info(genes):
    """ Retrieve start, end and strand data from gene name """
    matches = coll_unigenes.find({ 'o.g': {'$in': genes} })
    info = {}
    for m in matches:
        for gene_info in m.get('o', []):
            gene = gene_info['g']
            if gene not in genes:
                continue
            start, end, strand = gene_info['s']
            info[gene] = {
                'unigene': m['u'],
                'gene': gene, 
                'strand': strand, 
                'start': start,
                'end': end,
                'size': abs(int(end) - int(start)),
            }
    return info


def get_unigene_info(unigenes):
    functional_info = get_functional_info(unigenes)
    tax_preds = get_taxonomic_prediction(unigenes)
    info = {}
    info = { 
            unigene: {
                'tax_prediction': tax_preds.get(unigene, ''),
                **functional_info.get(unigene, {})
            } for unigene in unigenes }
    return info


def neighbor_analysis(unigenes):
    cluster_orfs = get_orfs(unigenes)
    print("\nObtaining neighbors\n")
    cluster_neighbors = { unigene: get_neighbors(orfs)
                  for unigene, orfs in cluster_orfs.items() }
    cluster_neighborhoods = {} 
    all_unigenes = set()

    print("\nComputing neighborhoods...\n")
    progress_idx = 0
    # First obtain most common contigs for each unigene in the cluster
    for unigene, orfs in cluster_neighbors.items():

        # PROGRESS
        printProgressBar(progress_idx, len(cluster_neighbors.keys()), 
                prefix = 'Progress:', length = 50)
        progress_idx += 1
        # PROGRESS

        neighborhoods = {}
        for v in orfs.values():
            central_strand, neigh_list = v.values()
            # Swap neighborhood order if central gene in negative strand
            # All central genes must be positively oriented
            if central_strand == '-':
                neigh_list.reverse()
            gene_info = get_gene_info(neigh_list)
            # Swap strands if central_strand is negative
            swap_strands(gene_info, central_strand)
            unigene_list = [g['unigene'] for g in gene_info.values()]
            # Order gene info by neighbor position
            gene_info_list = [gene_info.get(n, {'gene': n, 'strand': '+'})
                              for n in neigh_list]
            neighborhoods[tuple(unigene_list)] = gene_info_list
        # Most common neighborhood (by unigene assignation)
        n_keys = list(neighborhoods.keys())
        most_common = max(n_keys, key=n_keys.count)
        most_common_unigenes = set(most_common)
        most_common_genes = neighborhoods[most_common]

        all_unigenes.update(most_common_unigenes)
        cluster_neighborhoods[unigene] = most_common_genes

    # Time unigene info query
    t0 = time.time()
    unigene_info = get_unigene_info(list(all_unigenes))
    print(f'\n{round(time.time()-t0, 3)}s to query unigene info\n')

    cluster_neighborhood_info = []
    for unigene, neighborhood in cluster_neighborhoods.items():
        for i, neighbor_info in enumerate(neighborhood):
            n_info = {
                    'anchor': unigene, 
                    'pos': int(i - neighbor_range),
                    **neighbor_info,
                    **unigene_info.get(neighbor_info.get('unigene', ''), {}),
            }
            cluster_neighborhood_info.append(n_info)
    return cluster_neighborhood_info


def query_fam(query, n_range=10, cutoff=0):
    
    global client, db, coll_unigenes, coll_clusters, coll_e5, coll_taxa
    client,\
        db,\
        coll_unigenes,\
        coll_clusters,\
        coll_e5,\
        coll_taxa,\
        coll_seqs = mongo_connect()

    global ncbi
    ncbi = NCBITaxa()

    global max_gmgc_genes, percentage_cutoff, neighbor_range, max_orfs
    # maximun number of unigene allowed by gmgc cluster to be computed,
    # smaller the number smaller computing time
    max_gmgc_genes = int(400)
    # percentage of kegg/eggnog conservation in neighbor unigenes
    percentage_cutoff = float(cutoff)
    # number of neighbors genes to analyze around each unigene
    neighbor_range = n_range
    # max orfs per unigene to summarize neighborhood
    max_orfs = 5

    ### KEGG pathways
    global kegg_dict, eggNOG_DICT, egg_levels
    kegg_path = STATIC_PATH + "pickle/KEGG_DESCRIPTION.pickle"
    kegg_dict = get_pickle(kegg_path)
    egg_path = STATIC_PATH + "pickle/TAX_LEVELS.pickle"
    eggNOG_DICT = get_pickle(egg_path)
    # egg_levels = {}

    analysis_path = RESULTS_PATH + query + ".json"
    if path.exists(analysis_path):
        with open(analysis_path, "r") as infile:
            return load(infile)

    # Clean results tmp
    # rmtree(RESULTS_PATH)
    makedirs(RESULTS_PATH, exist_ok = True)

    ### Analysis
    analysis = {}

    member_list = write_newick(query, RESULTS_PATH)

    db_member_list = get_members(unigene_to_cl(query))
    print("\nSame members in mongodb and tree: "+
        str(sorted(member_list) == sorted(db_member_list)\
         and len(member_list) == len(db_member_list)))

    unigene_list = [clean_unigene(m) for m in member_list]

    print(f'\nNumber of members: {len(unigene_list)}\n')

    t0 = time.time()
    analysis = neighbor_analysis(unigene_list)
    print(f'\n{round(time.time()-t0, 3)}s to complete neighborhood analysis\n')

    with open(analysis_path, "w") as outfile:
        dump(analysis, outfile)

    return analysis
