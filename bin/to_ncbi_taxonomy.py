#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("lineage", help="File with one lineage string per line (e.g. dada2 output)")
parser.add_argument("names", help="The names.dmp file from the NCBI Taxonomy database")
args = parser.parse_args()

name_to_ncbi_id = {}
with open(args.names) as f:
    for line in f:
        (tax_id, name_txt, _, name_class, *_) = line.rstrip().replace('\t',' ').replace(' | ',' |').split(' |')
        if name_class == "scientific name":
            name_to_ncbi_id[name_txt] = tax_id

with open(args.lineage) as f:
    for line in f:
        taxa = line.rstrip().replace('; ',';').split(';')
        ncbi_id = 1
        species = ' '.join(taxa[-2:])
        if species in name_to_ncbi_id:
            ncbi_id = name_to_ncbi_id[species]
        else:
            for taxon in reversed(taxa[:-1]):
                if taxon in name_to_ncbi_id:
                    ncbi_id = name_to_ncbi_id[taxon]
                    break;
        print("{}".format(ncbi_id))
