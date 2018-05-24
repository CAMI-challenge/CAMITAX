#!/usr/bin/env python

import argparse
from collections import Counter

name_to_ncbi_id = {}
parent_id_map = {}
defined_ranks = {}
sequenced_taxa = set()

def getParent(ncbi_id):
    """Returns a node's parent node (or zero for the root node)"""
    parent_id = 0
    if ncbi_id in parent_id_map:
        parent_id = parent_id_map[ncbi_id]
    return parent_id if ncbi_id > 1 else 0

def getDepth(ncbi_id):
    """Returns the "real" depth of a taxon in the taxonomic tree"""
    while ncbi_id not in defined_ranks and ncbi_id > 1:
        ncbi_id = getParent(ncbi_id)
    return defined_ranks[ncbi_id] if ncbi_id > 1 else 0

def getLineage(ncbi_id):
    """Returns the defined full lineage (as NCBI IDs) of a taxon"""
    lineage = []
    while ncbi_id != 0:
        if ncbi_id in defined_ranks:
            lineage.insert(0, ncbi_id)
        ncbi_id = getParent(ncbi_id)
    return lineage

def getLCA(taxon_list):
    """Returns the lowest common ancestor (LCA) of all taxa in taxon_list"""
    if not taxon_list:
        return 1
    lineage = getLineage(taxon_list[0])
    for ncbi_id in taxon_list[1:]:
        while ncbi_id != 0:
            if ncbi_id in lineage:
                del lineage[lineage.index(ncbi_id)+1:]
                break
            ncbi_id = getParent(ncbi_id)
    return lineage[-1]

def getLowest(taxon_list):
    """Returns the lowest taxon in taxon_list, LCA of co-optimal solutions"""
    if not taxon_list:
        return 1
    taxon, depth = 1, 0
    for ncbi_id in taxon_list:
        t, d = ncbi_id, getDepth(ncbi_id)
        if d > depth:
            taxon, depth = t, d
        elif d == depth:
            taxon = getLCA([taxon, t])
    return taxon

def getIuLCA(taxon_list):
    """Returns the intersection union LCA of all taxa in taxon_list"""
    if not taxon_list:
        return 1
    ranks = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[]}
    for taxon in taxon_list:
        for ncbi_id in getLineage(taxon):
            ranks[defined_ranks[ncbi_id]].append(ncbi_id)
    iuLCA_taxon = 1
    for rank in range(7, 0, -1):
        rank_counter = Counter(ranks[rank])
        if rank_counter:
            iuLCA_taxon = getLowest([k for m in [min(rank_counter.values())] for k,v in rank_counter.items() if v >= 0.5*len(taxon_list)])
        if iuLCA_taxon > 1:
            break
    return iuLCA_taxon

def getRTLpath(taxon_list):
    """Returns the maximal root-to-leaf path of all taxa in taxon_list"""
    if not taxon_list:
        return 1
    ncbi_id_counter = Counter(taxon_list)
    if not ncbi_id_counter:
        ncbi_id_counter[1] += 1
    rtl_counter = Counter()
    for ncbi_id in ncbi_id_counter:
        leaf_taxon = ncbi_id
        node_taxon = ncbi_id
        while node_taxon != 0:
            node_count = ncbi_id_counter[node_taxon]
            rtl_counter[leaf_taxon] += node_count
            node_taxon = getParent(node_taxon)
    return getLCA([k for m in [max(rtl_counter.values())] for k,v in rtl_counter.items() if v == m])

def getLowNode(taxon_list):
    """Get the lowest node without siblings in the tree spanned by taxon_list"""
    if not taxon_list:
        return 1
    ranks = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[]}
    for taxon in taxon_list:
        for ncbi_id in getLineage(taxon):
            ranks[defined_ranks[ncbi_id]].append(ncbi_id)
    low_node = 1
    for rank in range(1, 8):
        rank_counter = Counter(ranks[rank])
        n = len(rank_counter)
        if n == 1:
            low_node = ranks[rank][0]
        elif n > 1:
            break
        # n == 0: incomplete lineage, skip gaps (e.g. for members of candidate phyla)
    return low_node

def getTaxID(name):
    """Match scientific name to NCBI Taxonomy ID"""
    taxa = name.rstrip().replace('; ',';').split(';')
    ncbi_id = 1
    species = ' '.join(taxa[-2:])
    if species in name_to_ncbi_id:
        ncbi_id = name_to_ncbi_id[species]
    else:
        for taxon in reversed(taxa[:-1]):
            if taxon in name_to_ncbi_id:
                ncbi_id = name_to_ncbi_id[taxon]
                break;
    return ncbi_id

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--names", help="NCBI Taxonomy file: names.dmp", required=True)
    parser.add_argument("--nodes", help="NCBI Taxonomy file: nodes.dmp", required=True)
    parser.add_argument("--mash", help="Mash output", required=True)
    parser.add_argument("--dada2", help="Dada2 lineage", required=True)
    parser.add_argument("--centrifuge", help="Centrifuge output", required=True)
    parser.add_argument("--kaiju", help="Kaiju output", required=True)
    parser.add_argument("--checkm", help="CheckM output", required=True)
    parser.add_argument("--known", help="taxIDs of sequenced genomes", required=True)
    parser.add_argument("--animax", help="best ANI hit", required=True)
    parser.add_argument("name", help="Genome name/prefix")
    args = parser.parse_args()

    with open(args.nodes) as f:
        ranks = {'superkingdom':1, 'phylum':2, 'class':3, 'order':4, 'family':5, 'genus':6, 'species':7}
        for line in f:
            (ncbi_id, parent_id, rank, *_) = line.split('|')
            ncbi_id, parent_id, rank = int(ncbi_id.strip()), int(parent_id.strip()), rank.strip()
            parent_id_map[ncbi_id] = parent_id
            if rank in ranks:
                defined_ranks[ncbi_id] = ranks[rank]

    with open(args.names) as f:
        for line in f:
            (ncbi_id, name, _, name_class, *_) = line.split('|')
            ncbi_id, name, name_class = int(ncbi_id.strip()), name.strip(), name_class.strip()
            if name_class == "scientific name":
                name_to_ncbi_id[name] = ncbi_id

    mash_taxon_list = []
    with open(args.mash) as f:
        for line in f:
            mash_taxon_list.append(int(line.strip()))
    mash_taxonomy = getLCA(mash_taxon_list)

    dada2_taxon_list = []
    with open(args.dada2) as f:
        for line in f:
            dada2_taxon_list.append(getTaxID(line.strip()))
    dada2_taxonomy = getIuLCA(dada2_taxon_list)

    centrifuge_taxon_list = []
    with open(args.centrifuge) as f:
        for line in f:
            centrifuge_taxon_list.append(int(line.strip()))
    centrifuge_taxonomy = getIuLCA(centrifuge_taxon_list)

    kaiju_taxon_list = []
    with open(args.kaiju) as f:
        for line in f:
            kaiju_taxon_list.append(int(line.strip()))
    kaiju_taxonomy = getIuLCA(kaiju_taxon_list)

    checkm_taxon_list = []
    with open(args.checkm) as f:
        next(f)
        for line in f:
            (_, name, _, _, _, completeness, contamination, strain_heterogeneity, genome_size, _, _, contigs, _, N50_contigs, _, _, _, _, gc_content, _, coding_density, _, predicted_genes, *_) = line.split('\t')
            sname = name[3:-(len(name)-name.index(" (UID"))]
            checkm_taxon_list.append(getTaxID(sname))
    checkm_taxonomy = getIuLCA(checkm_taxon_list)

    low_taxonomy = getLowNode([ mash_taxonomy, dada2_taxonomy,
                                centrifuge_taxonomy, kaiju_taxonomy,
                                checkm_taxonomy ])

    rtl_taxonomy = getRTLpath([ mash_taxonomy, dada2_taxonomy,
                                centrifuge_taxonomy, kaiju_taxonomy,
                                checkm_taxonomy ])

    with open(args.known) as f:
        for line in f:
            for taxon in getLineage(int(line.strip())):
                sequenced_taxa.add(taxon)

    with open(args.animax) as f:
        for line in f:
            ani_max = 100*float(line.strip()) # TODO Should be moved to parste_mash_output() eventually


    taxID = low_taxonomy

    novelty_ranks = {0:'NA', 1:'Known superkingdom', 2:'Known phylum', 3:'Known class', 4:'Known order', 5:'Known family', 6:'Known genus', 7:'Known species'}
    tax_lineage = getLineage(taxID)
    novelty_category = novelty_ranks[0]
    for node_taxon in tax_lineage[::-1]:
        if node_taxon in sequenced_taxa:
            novelty_category = novelty_ranks[defined_ranks[node_taxon]]
            break

    print("Genome\ttaxID\tNovelty_category\tANI\tCompleteness\tContamination\tContigs\tN50\tGC\tCoding_density\tPredicted_genes\tLowest\tRTLpath\tMash\tDada2\tCentrifuge\tKaiju\tCheckM\t")
    print("{}\t{}\t{}\t{:0.2f}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
        .format(args.name, taxID, novelty_category, ani_max, completeness, contamination,
                contigs, N50_contigs, gc_content, coding_density, predicted_genes,
                low_taxonomy, rtl_taxonomy, mash_taxonomy, dada2_taxonomy, centrifuge_taxonomy, kaiju_taxonomy, checkm_taxonomy))


if __name__ == "__main__":
    main()
