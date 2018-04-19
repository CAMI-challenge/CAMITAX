#!/usr/bin/env nextflow

params.i = 'input'
params.x = 'fasta'
genomes = Channel.fromPath("${params.i}/*.${params.x}")
db = Channel.fromPath( "$baseDir/db/", type: 'dir').first()

params.mash_db = 'rdp'
process mash {

    publishDir 'data', mode: 'copy'

    input:
    file db
    file genome from genomes

    output:
    file "${genome.baseName}.mash_10p.tsv"

    script:
    mash_index = "${db}/RefSeq_completeGenomes_20180410.msh"

    """
    mash dist ${mash_index} ${genome} | sort -gk3 | awk 'NR == 1 {t=\$3*1.1}; \$3 <= t' > ${genome.baseName}.mash_10p.tsv
    """
}

// process run_centrifuge {
//
//     publishDir 'data'
//
//     input:
//     file genes from prodigal_fna
//
//     output:
//     file "${genes.baseName}.centrifuge.txt"
//     file "${genes.baseName}.centrifuge.50"
//
//     """
//     ~/GitHub/centrifuge/centrifuge -f -x $baseDir/db/p_compressed ${genes} | ~/GitHub/centrifuge/centrifuge-kreport -x $baseDir/db/p_compressed > ${genes.baseName}.centrifuge.txt
//     cat ${genes.baseName}.centrifuge.txt | awk '+\$1 >= 50.0' | tail -n1 > ${genes.baseName}.centrifuge.50
//     """
// }
