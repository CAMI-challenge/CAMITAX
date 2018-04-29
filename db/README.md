# DADA2 taxonomic reference data

### Source: https://benjjneb.github.io/dada2/training.html

- Silva version 132: https://doi.org/10.5281/zenodo.1172782
- RDP trainset 16: https://doi.org/10.5281/zenodo.801827
- Silva DB might need filtering? Then host our own.

#  Mash sketch database

- We'll provide a more recent sketch of RefSeq release 87, insert instructions here:
  - In a nutshell: download RefSeq, rename FASTA files, sketch all.

# Kaiju index (also downloads NCBI taxonomy)

- Simply execute Kaiju's ``makeDB.sh -p``, which will download
  - the proGenomes database (``representatives.proteins.fasta.gz``), and
  - the NCBI Taxonomy (``names.dmp``, ``nodes.dmp``, and ``merged.dmp``).

# Centrifuge index

- Download the proGenomes database:
  - ``wget http://progenomes.embl.de/data/repGenomes/representatives.genes.fasta.gz``
- Replace taxon IDs found in ``merged.dmp`` by their updated IDs:
  - ``gunzip -c representatives.genes.fasta.gz | perl -lsne 'BEGIN{open(F,$m);while(<F>){@F=split(/[\|\s]+/);$h{$F[0]}=$F[1]}}if(/>(\d+)\.(\S+)/){print ">",defined($h{$1})?$h{$1}:$1,".",$2;}else{print}' -- -m=merged.dmp > centrifuge_db.fna``
- Create mapping file (tab-separated file mapping sequence IDs to taxon IDs):
  - ``perl -lne 'if(/>(\d+)\.(\S+)/){print $1,".",$2,"\t",$1}' < centrifuge_db.fna > centrifuge_db.conv``
- Build the index (creates ``proGenomes.[1-4].cf``)
  - ``centrifuge-build --conversion-table centrifuge_db.conv --taxonomy-tree nodes.dmp --name-table names.dmp centrifuge_db.fna proGenomes``
