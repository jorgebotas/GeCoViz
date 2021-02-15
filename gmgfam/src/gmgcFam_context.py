from pymongo import MongoClient
from ete3 import Tree
import gridfs

def mongo_connect():
    client = MongoClient('10.0.3.1', 27017, maxPoolSize=10)
    db = client.novel_fam
    gf = db.gene_families
    gmgcv1_neighs = db.gmgcv1_neighs
    # human_gut_neighs = db.neighs_human_gut
    # tara_mags_neighs = db.tara_mags_neighs
    # earth_mags_neighs = db.earth_mags_neighs
    # tara_euk_mags_neighs = db.tara_euk_MAGs
    return client,  gf, gmgcv1_neighs
           # human_gut_neighs, \
           # tara_mags_neighs, \
           # earth_mags_neighs, \
           # tara_euk_mags_neighs

def getDomains(query, client):
    """Retrieve pfam domains

    :query: unigene
    :db: MongoDB
    :returns: array of dicts with domain description

    """
    try:
        dms = client.gmgc_unigenes.pfam.find({
            'u' : query})[0]['pf']
        doms = []
        for d in dms:
            doms.append({
                     'class' : d['n'],
                     'start' : d['s'],
                     'end' : d['e'],
                     'shape' : 'rect'
                    })
    except:
        doms = False
    return doms


def formatContext(context, client):
    """Format context to fit new format

    :context: old format context
    :returns: new format context

    """
    newFormat = []
    for anchor, v in context.items():
        neighborhood = v['neighbourhood']
        for pos, neigh in neighborhood.items():
            gene = neigh['gene']
            geneName = neigh['code']
            if geneName == 'UNKNOWN': geneName = ''
            strand = neigh['strand']
            start = neigh['start']
            end = neigh['end']
            tax = neigh['tax_prediction']
            taxonomy = []
            for vals in tax.values():
                p = list(vals.values())[0]
                taxonomy.append({
                    'id' : p['id'],
                    'level' : p['rank'],
                    'description' : p['description']
                })
            ks = neigh['KEGG']
            kegg = []
            for id_, desc in ks.items():
                kegg.append({
                    'id' : id_,
                    'description' : desc['description']
                })
            eggs = neigh['eggNOG']
            eggnog = []
            for lev, val in eggs.items():
                try:
                    for id_, desc in val.items():
                        eggnog.append({
                            'id' : id_,
                            'level' : lev,
                            'description' : desc['description']
                        })
                except:
                    print(val)
            geneInfo = {
                'anchor' : anchor,
                'pos' : pos,
                'gene' : gene,
                'showName' : geneName,
                'strand' : strand,
                'start' : start,
                'end' : end,
                'kegg' : kegg,
                'eggnog' : eggnog,
                'taxonomic prediction' : taxonomy,
            }
            domains = getDomains(gene, client)
            if domains:
                geneInfo['domains'] = domains
            newFormat.append(geneInfo)
    return newFormat

def get_newick(query):
    client = mongo_connect()[0]
    treedb = client.trees
    fs = gridfs.GridFS(treedb)
    tfile = fs.find_one({"filename":query})
    treedata = fs.get(tfile._id)

    orig_newick = str(treedata.read(), "utf-8:")

    #Clean up newick
    t = Tree(orig_newick)
    for node in t:
        node.name = node.name.split('.')[1]
    newick = t.write()
    return newick

def get_context(query):
    client,  gf, gmgcv1_neighs = mongo_connect()
    gfam = query
    # gfam = gf.find({"gfn" : int(query)})[0]["gf"]
    context = gmgcv1_neighs.find({"gf" : int(gfam)})[0]['neigh']
    context = formatContext(context, client)
    return context
