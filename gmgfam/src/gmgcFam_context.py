from pymongo import MongoClient

def mongo_connect():
    client = MongoClient('10.0.3.1', 27017, maxPoolSize=10)
    db = client.novel_fam
    gf = db.gene_families
    gmgcv1_neighs = db.gmgcv1_neighs
    # human_gut_neighs = db.neighs_human_gut
    # tara_mags_neighs = db.tara_mags_neighs
    # earth_mags_neighs = db.earth_mags_neighs
    # tara_euk_mags_neighs = db.tara_euk_MAGs
    return gf, gmgcv1_neighs
           # human_gut_neighs, \
           # tara_mags_neighs, \
           # earth_mags_neighs, \
           # tara_euk_mags_neighs

def formatContext(context):
    """Format context to fit new format

    :context: old format context
    :returns: new format context

    """
    newFormat = []
    for anchor, v in context.items():
        neighborhood = v['neighbourhood']
        for pos, neigh in neighborhood.items():
            gene = neigh['gene']
            geneName = neigh['preferred_name']
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
                    id_, desc = val.items()
                    eggnog.append({
                        'id' : id_,
                        'level' : lev,
                        'description' : desc['description']
                    })
                except: pass
            geneInfo = {
                'anchor' : anchor,
                'pos' : pos,
                'gene' : gene,
                'showName' : geneName,
                'strand' : strand,
                'start' : start,
                'end' : end,
                'taxonomy' : taxonomy,
                'kegg' : kegg,
                'eggnog' : eggnog,
            }
            try:
                domains = neigh['domains']
                geneInfo['domains'] = domains
            except: pass
            newFormat.append(geneInfo)
    return newFormat

def get_context(query):
    gf, gmgcv1_neighs = mongo_connect()
    gfam = gf.find({"gfn" : int(query)})[0]["gf"]
    context = gmgcv1_neighs.find({"gf" : int(gfam)})[0]['neigh']
    context = formatContext(context)
    return context
