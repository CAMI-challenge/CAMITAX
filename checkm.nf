#!/usr/bin/env nextflow

params.i = 'input'
params.x = 'fasta'
params.db = "camitax/db"

Channel
    .fromPath("${params.i}/*.${params.x}")
    .map { file -> [file.baseName, file] }
    .into {
        input_genomes;
        checkm_genomes
    }

println("Input genomes:")
input_genomes.println { it[0] + " - " + it[1] }
db = Channel.fromPath( "${params.db}", type: 'dir').first()

process checkm {
    tag "${id}"

    publishDir "checkm/${id}"
    cpus = 8
    memory '40 GB'

    input:
    file db
    set val(id), file(genome) from checkm_genomes

    output:
    file "${id}.checkm.tsv"

    script:
    checkm_db = "${db}/checkm/"

    """
    # Set CheckM root data location
    echo ${checkm_db} | checkm data setRoot ${checkm_db}

    # Phylogenetic placement onto reduced reference tree
    checkm lineage_wf -t ${task.cpus} -x ${params.x} . checkm
    checkm qa -o 2 --tab_table checkm/lineage.ms checkm > ${id}.checkm.tsv
    """
}
