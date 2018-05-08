#!/usr/bin/env python

import argparse
from collections import Counter


parser = argparse.ArgumentParser()
parser.add_argument("taxon", help="File with one NCBI Taxonomy ID per line")
parser.add_argument("nodes", help="The nodes.dmp file from the NCBI Taxonomy database")
args = parser.parse_args()


parent_id_map = {}
with open(args.nodes) as f:
    for line in f:
        (ncbi_id, parent_id, *_) = line.split('|')
        ncbi_id, parent_id = int(ncbi_id.strip()), int(parent_id.strip())
        parent_id_map[ncbi_id] = parent_id
# Set the root node's parent ID to zero
parent_id_map[1] = 0


num_assignments = 0
ncbi_id_counter = Counter()
with open(args.taxon) as f:
    for line in f:
        num_assignments += 1
        ncbi_id = int(line.rstrip())
        ncbi_id_counter[ncbi_id] += 1


rtl_counter = Counter()
lca_counter = Counter()
for ncbi_id in ncbi_id_counter:
    leaf_taxon = ncbi_id
    leaf_count = ncbi_id_counter[ncbi_id]
    node_taxon = leaf_taxon
    while node_taxon != 0:
        node_count = ncbi_id_counter[node_taxon]
        rtl_counter[leaf_taxon] += node_count
        lca_counter[node_taxon] += leaf_count
        node_taxon = parent_id_map[node_taxon]
iulca_counter = { k:v for k, v in lca_counter.items() if v >= 0.5*num_assignments }


def getLineage(ncbi_id):
    lineage = []
    while ncbi_id != 0:
        lineage.insert(0, ncbi_id)
        ncbi_id = parent_id_map[ncbi_id]
    return lineage

def getLCA(taxon_list):
    lineage = getLineage(taxon_list.pop())
    for ncbi_id in taxon_list:
        while ncbi_id != 0:
            if ncbi_id in lineage:
                del lineage[lineage.index(ncbi_id)+1:]
                break
            ncbi_id = parent_id_map[ncbi_id]
    return lineage.pop()

def getLowest(taxon_list):
    depth = 0
    taxon = 0
    for ncbi_id in taxon_list:
        d = 0
        t = ncbi_id
        while ncbi_id != 0:
            d += 1
            ncbi_id = parent_id_map[ncbi_id]
        if d > depth:
            depth = d
            taxon = t
        elif d == depth:
            taxon = getLCA([taxon, t])
    return taxon

lca_support = lca_counter[max(lca_counter, key=lca_counter.get)]
lca_candidates = [ k for k, v in lca_counter.items() if v == lca_support ]
lca_winner = getLCA(list(ncbi_id_counter.keys()))

iulca_support = iulca_counter[min(iulca_counter, key=iulca_counter.get)]
iulca_candidates = [ k for k, v in iulca_counter.items() if v == iulca_support ]
iulca_winner = getLowest(iulca_candidates.copy())

rtl_support = rtl_counter[max(rtl_counter, key=rtl_counter.get)]
rtl_candidates = [ k for k, v in rtl_counter.items() if v == rtl_support ]
rtl_winner = getLCA(rtl_candidates.copy())

print("#\tTaxon\tSupport\tCandidates")
print("LCA\t{}\t{}\t{}".format(lca_winner, lca_support, lca_candidates))
print("iuLCA\t{}\t{}\t{}".format(iulca_winner, iulca_support, iulca_candidates))
print("RTLpath\t{}\t{}\t{}".format(rtl_winner, rtl_support, rtl_candidates))
