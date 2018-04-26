#!/usr/bin/env nextflow

params.i = 'input'
params.x = 'fasta'

Channel
    .fromPath("${params.i}/*.${params.x}")
    .into {
        input_genomes;
        barrnap_genomes;
        bedtools_genomes;
        prodigal_genomes;
        mash_genomes;
        checkm_genomes
    }

input_genomes.collect().println { "Input genomes: " + it }
db_folder = Channel.fromPath( "${baseDir}/db/", type: 'dir').first()

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

checkm_in = prodigal_faa.collect()
process checkm {
    publishDir 'data', mode: 'copy'
    maxForks 1

    input:
    file db_folder
    file checkm_in

    output:
    file "checkm.qa"

    """
    echo ${db_folder} | checkm data setRoot ${db_folder}
    checkm lineage_wf --reduced_tree --genes -x faa . checkm_out
    checkm qa -o 2 --tab_table checkm_out/lineage.ms checkm_out > checkm.qa
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
