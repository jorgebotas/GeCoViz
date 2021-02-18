from pymongo import MongoClient

def mongoConnect():
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
