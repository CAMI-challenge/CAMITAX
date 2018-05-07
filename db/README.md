## CheckM database

- **Source:** https://ecogenomics.github.io/CheckM/
  - https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm
  - https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz


## DADA2 taxonomic reference data [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.801828.svg)](https://doi.org/10.5281/zenodo.801828)


- **Source:** https://benjjneb.github.io/dada2/training.html
  - Silva version 132: https://doi.org/10.5281/zenodo.1172782
  - RDP trainset 16: https://doi.org/10.5281/zenodo.801827


## Kaiju index (and NCBI taxonomy)

```
# Kaiju's makeDB.sh downloads the proGenomes and NCBI Taxonomy databases, and builds the index:
makeDB.sh -p
```

## Centrifuge index

```
# Download the proGenomes database (genes):
wget http://progenomes.embl.de/data/repGenomes/representatives.genes.fasta.gz

# Replace taxon IDs found in merged.dmp by their updated IDs:
gunzip -c representatives.genes.fasta.gz | perl -lsne 'BEGIN{open(F,$m);while(<F>){@F=split(/[\|\s]+/);$h{$F[0]}=$F[1]}}if(/>(\d+)\.(\S+)/){print ">",defined($h{$1})?$h{$1}:$1,".",$2;}else{print}' -- -m=merged.dmp > centrifuge_db.fna

# Generate mapping file (tab-separated file mapping sequence IDs to taxon IDs):
perl -lne 'if(/>(\d+)\.(\S+)/){print $1,".",$2,"\t",$1}' < centrifuge_db.fna > centrifuge_db.conv

# Build the index (creates proGenomes.[1-4].cf):
centrifuge-build --conversion-table centrifuge_db.conv --taxonomy-tree nodes.dmp --name-table names.dmp centrifuge_db.fna proGenomes
```


##  Mash sketch database

```
# Download (or update) all bacterial and archaeal genome sequences from RefSeq
rsync --exclude='*_from_genomic.*' --include='archaea/**/latest**/*_genomic.*' --include='bacteria/**/latest**/*_genomic.*' --include='*/' --exclude='*' --recursive --copy-links --prune-empty-dirs --times --verbose --delete-after rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/ refseq

# Use https://github.com/ondovb/refseqCollate to collate them
collateGenomes.sh refseq mash

# Create the Mash index
cd mash.genomes && mash sketch -l <(find -name "*.fasta" | cut -c3-) -o RefSeq
```
