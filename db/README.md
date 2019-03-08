##  Mash sketch database

```
# Follow the Genomes Download FAQ (minus the 'complete' part): https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#allcomplete
wget -O assembly_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
awk -F "\t" '$11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
xargs -n 1 -P 32 wget -q < ftpfilepaths

# Double-check that all genomes were downloaded, please manually download missing file(s) – if any:
diff <(sed 's/.*\///' < ftpfilepaths | sort) <(ls -1 | grep "genomic.fna.gz" | sort)

# Prefix the filename with the NCBI Taxonomy ID (field 6 in assembly_summary.txt):
mkdir mash_genomes && cut -f6,20 assembly_summary.txt | sed 's/ftp.*\///' | while read i; do cp $(echo $i | cut -d' ' -f2)\_genomic.fna.gz mash_genomes/$(echo $i | tr ' ' '_'); done

# Create the Mash index
cd mash_genomes && mash sketch -p 32 -o ../RefSeq * && cd ..
mash info -t RefSeq.msh | tail -n +2 | cut -f3 | cut -f1 -d'_' | sort | uniq > RefSeq.ids
```

## DADA2 reference databases

```
# Download RDP trainset 16: https://doi.org/10.5281/zenodo.801827
wget https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz
wget https://zenodo.org/record/801828/files/rdp_species_assignment_16.fa.gz

# Download SILVA version 132: https://doi.org/10.5281/zenodo.1172782
wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz
wget https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz

# Source: https://benjjneb.github.io/dada2/training.html
```

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

## CheckM database

```
# Download and unpack CheckM database
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar xfv checkm_data_2015_01_16.tar.gz

# Source: https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm
```
