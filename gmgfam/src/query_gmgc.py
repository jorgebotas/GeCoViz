from collections import defaultdict
from ete3 import Tree, NCBITaxa
import gridfs
from os import makedirs
import pickle
from pymongo import MongoClient
from shutil import rmtree

from django.conf import settings

RESULTS_PATH = settings.BASE_DIR + '/gmgfam/results.tmp/'
STATIC_PATH = settings.BASE_DIR + '/static/gecoviz/'


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

    return [client, db, coll_unigene, coll_clusters, coll_e5, coll_taxa]


def get_pickle(filepath):
    """generate kegg pathway dictionary containning kegg descriptions"""
    # read dictionary from pickle file
    with open(filepath, 'rb') as pickle_in:
        pdict = pickle.load(pickle_in)
    return pdict


def clean_unigene(gmgc):
        '''in case GMGC nomeclature was used in the input example,
         GMGC.100_000_123.UNKNOWN words after/before . are removed'''
        gmgc_clean = gmgc.split(".")[1]
        return gmgc_clean


def cl_to_unigene(cl):
    u = coll_clusters.find_one({'cl' : cl}).get('u', cl)
    return u


def unigene_to_cl(unigene):
    cl = coll_clusters.find_one({'u' : unigene}).get('cl', '')
    return cl


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
                'prefered_name': cl.get('p_n', ''),
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
    with open(results_dir + query_cluster + "_newick.txt", "w") as outputfile:
        outputfile.write(newick)
    return


def swap_strand(s, reference_s):
    if reference_s == "+":
        return s
    else:
        if s == "+":
            return "-"
        elif s == "-":
            return "+"
    return "NA"


def get_orfs(unigenes):
    """ Retrieve orf genomic information
    Returns dict(unigene: dict(central_gene: info))
    """
    cluster = coll_unigenes.find({ 'u': {'$in': unigenes} })
    orfs = defaultdict(dict)
    for orf in cluster:
        unigene = orf['u']
        orf = orf['o']
        for a in orf:
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
    all_genes = set()
    for m in matches:
        gene_info = m.get('o', {})
        gene = gene_info['g']
        all_genes.add(gene)
        start, end, strand = gene_info['s']
        info[gene] = {
            'gene': gene, 
            'strand': strand, 
            'start': start,
            'end': end,
            'size': abs(int(end) - int(start)),
        }
    all_genes = list(all_genes)
    functional_info = get_functional_info(all_genes)
    prefered_names = get_preferred_name(all_genes)
    tax_preds = get_taxonomic_prediction(unigenes)
    for gene, gene_info in info.items():
        unigene_info = {
                'tax_prediction': tax_preds.get(gene, ''),
                'domains': tax_preds.get(gene, ''),
        }
        fn_info = functional_info.get(gene, {})
        info[gene] = {**gene_info, **unigene_info, **fn_info}

    return info


def neighbor_analysis(unigenes):
    cluster_orfs = get_orfs(unigene)
    cluster_neighbors = { unigene: get_neighbors(orfs)
                  for unigene, orfs in cluster_orfs.items() }
    all_genes = []

    # First obtain most common contigs for each unigene in the cluster
    for unigene, orfs in cluster_neighbors.items():
        neighborhoods = []
        for v in orfs.values():
            central_strand, neigh_list = v.values()
            # Swap neighborhood order if central gene in negative strand
            # All central genes must be positively oriented
            if central_strand == '-':
                neigh_list.reverse()
            neighborhoods.append(neigh_list)
        # Most common neighborhood
        most_common = max(neighborhoods, key=neighborhoods.count())
        cluster_neighborhoods[unigene] = most_common
        all_genes += most_common

    # Time gene info query
    t0 = time.time()
    gene_info = get_gene_info(all_genes)
    print(f'\n{time.time()-t0}s to query gene info\n')

    cluster_neighborhood_info = {}
    for unigene, neighborhood in cluster_neighborhoods.items():
        for i, neighbor in enumerate(neighborhood):
            n_info = {
                    'anchor': unigene, 
                    'pos': int(i - neighbor_range),
                    'gene': neighbor, 
                    **gene_info.get(neighbor, {}),
            }
            cluster_neighborhood_info.append(n_info)
    return cluster_neighborhood_info


def query_fam(query, n_range=2, cutoff=0):
    
    global client, db, coll_unigenes, coll_clusters, coll_e5
    client,\
        db,\
        coll_unigenes,\
        coll_clusters,\
        coll_e5,\
        coll_taxa = mongo_connect()

    global ncbi
    ncbi = NCBITaxa()

    global max_gmgc_genes, percentage_cutoff, neighbor_range
    # maximun number of unigene allowed by gmgc cluster to be computed,
    # smaller the number smaller computing time
    max_gmgc_genes = int(400)
    # percentage of kegg/eggnog conservation in neighbor unigenes
    percentage_cutoff = float(cutoff)
    # number of neighbors genes to analyze around each unigene
    neighbor_range = n_range

    ### KEGG pathways
    global kegg_dict, eggNOG_DICT, egg_levels
    kegg_path = STATIC_PATH + "pickle/KEGG_DESCRIPTION.pickle"
    kegg_dict = get_pickle(kegg_path)
    egg_path = STATIC_PATH + "pickle/TAX_LEVELS.pickle"
    eggNOG_DICT = get_pickle(egg_path)
    # egg_levels = {}

    # Clean results tmp
    rmtree(RESULTS_PATH)
    makedirs(RESULTS_PATH, exist_ok = True)

    ### Analysis
    analysis = {}

    write_newick(query, RESULTS_PATH)
    member_list = get_members(unigene_to_cl(query))
    unigene_list = [clean_unigene(cl_to_unigene(m)) for m in member_list]

    analysis = neighbor_analysis(unigene_list)

