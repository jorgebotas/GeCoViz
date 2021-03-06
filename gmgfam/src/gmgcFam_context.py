from django.conf import settings
from ete3 import Tree
import gridfs

from .mongoConnect import mongoConnect

RESULTS_PATH = settings.BASE_DIR + '/gmgfam/tmp/'

def getDomains(query, client):
    """Retrieve pfam domains

    :query: unigene
    :db: MongoDB
    :returns: array of dicts with domain description

    """
    try:
        dms = client.gmgc_unigenes.pfam.find({'u' : query})[0]['pf']
        doms = []
        for d in dms:
            doms.append({
                     'id' : d['n'],
                     'start' : d['s'],
                     'end' : d['e'],
                     'shape' : 'rect',
                     'description' : ''
                    })
    except:
        doms = False
    return doms

def getSequence(query, client):
    """Retrieve query sequence

    :query: unigene
    :client: MongoDB client
    :returns: sequence
    """
    try:
        seq = client.gmgc_unigenes.sequences.find({'u' : query})[0]['sq']
    except:
        seq = False
    return seq

def orderTaxonomy(taxonomy):
    """Orders taxonomy array of dicts
    by taxonomic rank

    :taxonomy: array of dicts. {id:str, description...}
    :returns: sorted array
    """
    taxRanks = [
        'no rank',
        # 'domain',
        'superkingdom',
        # 'kingdom',
        'phylum',
        'class',
        # 'superorder',
        'order',
        # 'superfamily',
        'family',
        'genus',
        'species',
    ]
    def orderKey(t):
        return taxRanks.index(t['level'])
    return sorted(taxonomy, key=orderKey)

def formatContext(context, client):
    """Format context to fit new format

    :context: old format context
    :returns: new format context

    """
    newFormat = []
    membersTaxonomy = {}
    for anchor, v in context.items():
        neighborhood = v['neighbourhood']
        membersTaxonomy[anchor] = []
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
            taxonomy = orderTaxonomy(taxonomy)
            if len(taxonomy) > len(membersTaxonomy[anchor]):
                membersTaxonomy[anchor] = taxonomy
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
                'gene name' : geneName,
                'strand' : strand,
                'start' : start,
                'end' : end,
                'kegg' : kegg,
                'eggnog' : eggnog,
                'taxonomy' : taxonomy,
            }
            seq = getSequence(gene, client)
            if seq:
                geneInfo['length'] = len(seq)
            domains = getDomains(gene, client)
            if domains:
                geneInfo['pfam'] = domains
            newFormat.append(geneInfo)
    return newFormat, membersTaxonomy

def get_newick(query, membersTaxonomy):
    client = mongoConnect()[0]
    treedb = client.trees
    fs = gridfs.GridFS(treedb)
    tfile = fs.find_one({"filename":query})
    treedata = fs.get(tfile._id)
    orig_newick = str(treedata.read(), "utf-8:")
    #Clean up newick
    t = Tree(orig_newick)
    for node in t:
        name = node.name.split('.')[1]
        showName = name
        taxonomy = membersTaxonomy[name]
        if len(taxonomy) > 0:
            for tax in taxonomy:
                name += '.' + tax['description'].strip().replace('.', '')
        if len(taxonomy) > 1:
            showName = taxonomy[-1]['description'].strip().replace('.', '')
        name = showName + '.' + name
        node.name = str(name)
    newick = t.write()
    return newick

def get_context(query):
    client, gf, gmgcv1_neighs = mongoConnect()
    if len(query.strip()) == 11:
        gfam = query
    else:
        gfam = gf.find({"gfn" : int(query)})[0]["gf"]
        gfam = str(gfam).zfill(9)
        gfam = '{}_{}_{}'.format(gfam[:3], gfam[3:6], gfam[6:])
    context = gmgcv1_neighs.find({"gf" : int(gfam)})[0]['neigh']
    context, membersTaxonomy = formatContext(context, client)
    newick = get_newick(str(gfam), membersTaxonomy)
    with open(RESULTS_PATH + str(query) + '_tree.nwx', "w") as handle:
        handle.write(newick)
    return context
