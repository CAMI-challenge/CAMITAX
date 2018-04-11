#!/usr/bin/env nextflow

params.i = 'input'
params.x = 'fasta'
genomes = Channel.fromPath("${params.i}/*.${params.x}")

process run_prodigal {

	container = 'quay.io/biocontainers/prodigal:2.6.3--0'

	publishDir 'data'

	input:
	file genome from genomes

	output:
	file "${genome.baseName}.prodigal.faa"
    file "${genome.baseName}.prodigal.fna" into prodigal_fna
    file "${genome.baseName}.prodigal.gff"

	"""
    prodigal -a ${genome.baseName}.prodigal.faa -d ${genome.baseName}.prodigal.fna -f gff -i ${genome} -o ${genome.baseName}.prodigal.gff
	"""
}

process run_centrifuge {

	publishDir 'data'

	input:
	file genes from prodigal_fna

	output:
    file "${genes.baseName}.centrifuge.txt"
    file "${genes.baseName}.centrifuge.50"

	"""
    ~/GitHub/centrifuge/centrifuge -f -x $baseDir/db/p_compressed ${genes} | ~/GitHub/centrifuge/centrifuge-kreport -x $baseDir/db/p_compressed > ${genes.baseName}.centrifuge.txt
    cat ${genes.baseName}.centrifuge.txt | awk '+\$1 >= 50.0' | tail -n1 > ${genes.baseName}.centrifuge.50
	"""
}
