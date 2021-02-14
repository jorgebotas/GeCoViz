from django.conf import settings
from csv import DictReader
from ete3 import Tree
from glob import glob
import json
import os.path
import re

from .mongoClient import mongoConnect

def getPickle(filePath):
    """
    Return dict contained in pickle file

    :filePath: path to pickle file
    :returns: dictionary
    """
    with open(filePath, 'rb') as pickle_in:
        pdict = pickle.load(pickle_in)
    return pdict

def get_eggDescription(eggnog, client):
    """
    Get eggnog OG description

    :eggnog: eggNOG id (OG)
    :db: MongoDB
    :returns: string with description,
              empty string if not found
    """
    try:
        description = client.gmgc_unigenes.eggnog_v5.find({"e" : eggnog})[0]["d"]
    except:
        description = ""
    return description

def get_eggNOG_tree(query, tree_file, db):
    """
    Retrieve eggNOG OG's tree from MongoDB
    Store Newick in tree_file

    :query: eggNOG id (OG)
    :tree_file: Newick output file
    :db: MongoDB
    """
    t = db.trees.find({"e" : query})[0]["tree"]
    tree = Tree(t)
    try:
        taxDict = getPickle(STATIC_PATH + "pickle/TAX_LEVELS.pickle")
    except:
        taxDict = False
    for node in tree:
        name = str(node.name).split(".")
        try: tax = taxDict[name[0]]
        except:  tax = name[0]
        tax = tax.replace(".", "").strip().split(" ")
        tax_s = "_".join(tax[:2])
        tax_l = "_".join(tax)
        node.name = ".".join([tax_s] + name[:1] + [tax_l] + name[1:])
    with open(tree_file, "w") as handle:
        handle.write(tree.write())
    return

def getMembers(query, db):
    """
    Retrieve members from specified OG

    :query: eggNOG id (OG)
    :db: MongoDB
    :returns: list with members
    """
    og = db.members.find({"e" : query})
    members = set() # members set
    for member in og:
        orf = member["s"]
        members.add(orf)
    return list(members)

def getNeighbors(members, nNeigh, db):
    """
    Retrieve neighborhood data for members of eggNOG OG

    :members: list of members
    :nNeigh: number of neigbors up and downstream
    :db: MongoDB
    :returns: dict where key is member and value a list of ordered neighbors
    """
    neighDict = {}
    for m in members:
        ordered_genes = []
        try:
            contig, orf, *_ = m.split("_")
            contig, orf = str(contig), str(orf)
        except: continue
        if "gene" in m:
            seqs = db.contig_clusters.find({"contig" : contig})[0]["o"]
        for pos in range(-nNeigh, nNeigh + 1):
            try:
                # Split letters and numbers
                prefix, number = re.search('([^0-9]+)(\d+)$', orf).groups()
                gene = prefix + str(int(number) + int(pos))
            except:
                gene = str(int(orf) + int(pos))
            if "gene" in m:
                if gene in seqs:
                    gene = "_".join([contig, gene])
                else:
                    # Not present in contig
                    gene = "NIcontig"
            else:
                if len(orf) != len(gene):
                    gene = gene.zfill(len(orf))
                gene = contig + str(gene)
            ordered_genes.append(gene)
        neighDict[m] = ordered_genes
    return neighDict

def getGeneData(gene, client, db, taxDict, keggDict):
    """
    Retrieve gene data

    :gene: gene name
    :returns: list []
    """
    # geneDesc = list(db.gene_description.find({"o" : gene}))
    geneDesc = []
    if len(geneDesc) > 0: geneDesc = geneDesc[0]["d"]
    else: geneDesc = ""
    geneInfo_fromContigs = list(db.contigs.find({"o" : gene}))
    if len(geneInfo_fromContigs) > 0:
        strand = geneInfo_fromContigs[0]["str"]
        start = geneInfo_fromContigs[0]["s"]
        end = geneInfo_fromContigs[0]["e"]
        orfID_forAnnotation = geneInfo_fromContigs[0]["o_i"]
        geneInfo_fromAnnotation = list(db.annotation.find({"o" :
                                        orfID_forAnnotation}))
        if len(geneInfo_fromAnnotation) > 0:
            geneName = geneInfo_fromAnnotation[0]["g_n"]
            geneDesc = geneInfo_fromAnnotation[0]["d"]
            kegg = geneInfo_fromAnnotation[0]["Kegg"]
            keggJSON = []
            if kegg != "":
                keggList = kegg.split(",")
                for k in keggList:
                    try:
                        desc = keggDict[k]
                    except:
                        desc = ""
                    keggJSON.append({
                        'id' : k,
                        'description' : desc
                    })
            eggnog = geneInfo_fromAnnotation[0]["enog"]
            eggJSON = []
            if eggnog != "":
                eggList = eggnog.split(",")
                for e in eggList:
                    e = e.split('@')
                    if taxDict: level = get_taxLevel(e[1], taxDict)
                    else: level = e[1]
                    eggJSON.append({
                        'id' : e[0],
                        'level' : level,
                        'description' : get_eggDescription(e[0], client)
                    })
        else:
            geneName, keggJSON, eggJSON = "", "", ""
        taxonomy, gene = gene.split('.');
        geneInfo = [gene, geneName, geneDesc,
                    strand, start, end,
                    taxonomy,
                    keggJSON, eggJSON]
    else:
        # Gene info not found
        taxonomy, gene = gene.split('.');
        geneInfo = [gene, "", "", "", "", "", taxonomy,  "", ""]
    return geneInfo

def getNeighData(neighDict, client, db):
    """
    Retrieves data for neighbors in neighDict gene neighborhoods

    :neighDict: dict where key is member and value a list of ordered neighbors
    :db: MongoDB
    :returns:
    """
    neighData = {}
    try:
        taxDict = getPickle(STATIC_PATH + "pickle/TAX_LEVELS.pickle")
    except:
        taxDict = False
    try:
        keggDict = getPickle(STATIC_PATH + "pickle/KEGG_DESCRIPTION.pickle")
    except:
        keggDict = False
    for member, neighs in neighDict.items():
        neighList = []
        for n in neighs:
            geneData = getGeneData(n, client, db, taxDict, keggDict)
            neighList.append(geneData)
        neighData[member] = neighList

    return neighData

def tsvToJson(filePath):
    """
    Retrieve tsv file as json format

    :filePath: path to .tsv
    :returns: json format of tsv file content

    """
    jsonData = []
    with open(filePath, 'r') as handle:
        tsvData = DictReader(handle, delimiter="\t")
        for line in tsvData:
            jsonData.append(line)
    return jsonData

def launch_analysis(query, nNeigh, tmpDir=False):
    """
    Generate tsv with gene neighborhood data
    and retrieve eggNOG OG tree (.nwx).
    Store in results directory

    :query: eggNOG id (OG)
    :nNeigh: number of neighborhood genes at each side
    """
    if not tmpDir:
        tmpDir = settings.BASE_DIR +  "/progenomes/tmp/"
    cacheLimit = 50 # Max number of files in tmpDir
    query = query.split('@')[0].strip()
    neighData_file = tmpDir + "{}_{}_neighData.tsv".format(query, nNeigh)
    fileList = glob(tmpDir + '*.tsv') + glob(tmpDir + '*.nwx')
    # Sort files by modification date
    fileList.sort(key=os.path.getmtime, reverse=True)
    # No need to recompute if file already in "cache"
    # if os.path.isfile(neighData_file): return tsvToJson(neighData_file)
    # Remove files in tmpDir if too many
    while len(fileList) > cacheLimit:
        os.remove(tmpDir + fileList.pop())
    # ANALYSIS
    client, db = mongoConnect()
    # Get tree
    tree_file = tmpDir + "{}_tree.nwx".format(query)
    # try:
    get_eggNOG_tree(query, tree_file, db)
    # except:
        # print("\nNo tree found for {}\n".format(query))
    # Get neighborhood data
    members = getMembers(query, db)
    neighDict = getNeighbors(members, nNeigh, db)
    neighData = getNeighData(neighDict, client, db)
    headers = ["gene", "anchor", "pos",
               "showName",  "description",
               "strand", "start", "end",
               "taxonomy",
              "kegg", "eggnog"]
    with open(neighData_file, "w") as handle:
        handle.write("\t".join(headers) + "\n")
        for anchor, neighborhood in neighData.items():
            for i, neigh in enumerate(neighborhood):
                pos = str(i - nNeigh)
                n, geneName, geneDesc,\
                    strand, start, end,\
                    taxonomy,\
                    keggJSON, eggJSON = neigh
                # Only keep gene identifier
                anchor = anchor.split('.')[-1]
                line = [n, anchor, pos,
                        geneName, geneDesc,
                        strand, start, end,
                        str(taxonomy),
                        str(keggJSON), str(eggJSON)]
                handle.write("\t".join(line) + "\n")
    return tsvToJson(neighData_file)

# with open('./tmp/43PAE_2_neighData.json', 'w') as handle:
    # jsonDump = json.dumps(launch_analysis("43PAE", 2),indent=2)
    # handle.write(jsonDump)
