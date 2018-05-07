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
        (ncbi_id, parent_id, rank, *_) = line.rstrip().replace('\t',' ').replace(' | ',' |').split(' |')
        parent_id_map[int(ncbi_id)] = int(parent_id)
# Set the root node's parent ID to zero
parent_id_map[1] = 0

# TODO Store each line in a set to compute union/intersection
num_assignments = 0
ncbi_id_counter = Counter()
with open(args.taxon) as f:
    for line in f:
        ncbi_id = line.rstrip()
        ncbi_id_counter[int(ncbi_id)] += 1
        num_assignments += 1

rtl_counter = Counter()
lca_counter = Counter()
for ncbi_id in ncbi_id_counter:
    id_taxon = ncbi_id
    id_count = ncbi_id_counter[id_taxon]
    while ncbi_id != 0:
        if ncbi_id in ncbi_id_counter:
            rtl_counter[id_taxon] += ncbi_id_counter[ncbi_id]
        lca_counter[ncbi_id] += id_count
        ncbi_id = parent_id_map[ncbi_id]
iulca_counter = { k:v for k, v in lca_counter.items() if v >= 0.5*num_assignments }

# TODO Take care of co-optimal solutions!!!
# TODO LCA should be deepest taxon in the list of co-optimal results
lca_winner = max(lca_counter, key=lca_counter.get)
# TODO iuLCA and mRTLp should be the LCA of co-optimal solutions
iulca_winner = min(iulca_counter, key=iulca_counter.get)
rtl_winner = max(rtl_counter, key=rtl_counter.get)
#max_winner = max(ncbi_id_counter, key=ncbi_id_counter.get)

print("#Metric\tValue\tSupport")
print("LCA\t{}\t{}".format(lca_winner, lca_counter[lca_winner]))
print("iuLCA\t{}\t{}".format(iulca_winner, iulca_counter[iulca_winner]))
print("mRTLp\t{}\t{}".format(rtl_winner, rtl_counter[rtl_winner]))
#print("MAX\t{}\t{}".format(max_winner, ncbi_id_counter[max_winner]))
