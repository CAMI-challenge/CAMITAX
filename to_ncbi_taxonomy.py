#!/usr/bin/env python3

import gzip
name_to_ncbi_id = {}
ncbi_id_to_lineage = {}

with gzip.open("/Users/abremges/GitHub/CAMITAX/db/taxonomy/rankedlineage.dmp.gz") as f:
	for line in f:
		line = ''.join(line.decode('utf-8').rstrip().split('\t'))
		taxa = line.split('|')
		del taxa[8] # Remove (empty) kingdom entry to match CAMI spec
		name_to_ncbi_id[taxa[1]] = taxa[0]
		ncbi_id_to_lineage[taxa[0]] = ';'.join((taxa[8:1:-1]))

with open("/Users/abremges/GitHub/CAMITAX/data/s1.16S_rRNA.dada2.txt") as f:
	with open("out.txt", 'w') as g:
		for line in f:
			line = line.rstrip()
			taxa = line.split(';')
			ncbi_id = 1
			if ' '.join(taxa[-2:]) in name_to_ncbi_id:
				ncbi_id = name_to_ncbi_id[' '.join(taxa[-2:])]
			else:
				for taxon in reversed(taxa[:-1]):
					if taxon in name_to_ncbi_id:
						ncbi_id = name_to_ncbi_id[taxon]
						break;
			g.write("{}\t{}\n".format(ncbi_id, ncbi_id_to_lineage[ncbi_id]))
