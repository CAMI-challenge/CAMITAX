#!/usr/bin/env nextflow

params.i = 'input'
params.x = 'fasta'

Channel
    .fromPath("${params.i}/*.${params.x}")
    .into {
        input_genomes;
        prodigal_genomes;
        diamond_genomes
    }

input_genomes.collect().println { "Input genomes: " + it }
db = Channel.fromPath( "${baseDir}/db/", type: 'dir').first()

process prodigal {
    publishDir 'data', mode: 'copy'

    input:
    file genome from prodigal_genomes

    output:
    file "${genome.baseName}.prodigal.faa" into prodigal_faa
    file "${genome.baseName}.prodigal.fna" into prodigal_fna

    """
    prodigal -a ${genome.baseName}.prodigal.faa -d ${genome.baseName}.prodigal.fna -f gff -i ${genome} -o ${genome.baseName}.prodigal.gff
    """
}


process diamond {
    publishDir 'data', mode: 'copy'
    maxForks 1

    input:
    file db
    file genome from diamond_genomes
    file "${genome.baseName}.prodigal.faa" from prodigal_faa

    output:
    file "${genome.baseName}.diamond.out"

    script:
    diamond_db = "${db}/diamond/swissprot_20180410.dmnd"

    """
    diamond blastp --sensitive --algo 1 -f 102 -p 1 -d ${diamond_db} -q ${genome.baseName}.prodigal.faa > ${genome.baseName}.diamond.out
    """
}


//

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
