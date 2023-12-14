#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
camitax_version = 'v0.7.0'
println """
 ▄████▄   ▄▄▄       ███▄ ▄███▓ ██▓▄▄▄█████▓ ▄▄▄      ▒██   ██▒
▒██▀ ▀█  ▒████▄    ▓██▒▀█▀ ██▒▓██▒▓  ██▒ ▓▒▒████▄    ▒▒ █ █ ▒░
▒▓█    ▄ ▒██  ▀█▄  ▓██    ▓██░▒██▒▒ ▓██░ ▒░▒██  ▀█▄  ░░  █   ░
▒▓▓▄ ▄██▒░██▄▄▄▄██ ▒██    ▒██ ░██░░ ▓██▓ ░ ░██▄▄▄▄██  ░ █ █ ▒
▒ ▓███▀ ░ ▓█   ▓██▒▒██▒   ░██▒░██░  ▒██▒ ░  ▓█   ▓██▒▒██▒ ▒██▒
░ ░▒ ▒  ░ ▒▒   ▓▒█░░ ▒░   ░  ░░▓    ▒ ░░    ▒▒   ▓▒█░▒▒ ░ ░▓ ░
  ░  ▒     ▒   ▒▒ ░░  ░      ░ ▒ ░    ░      ▒   ▒▒ ░░░   ░▒ ░
░          ░   ▒   ░      ░    ▒ ░  ░        ░   ▒    ░    ░
░ ░            ░  ░       ░    ░                 ░  ░ ░    ░
░                                    ${camitax_version?:''}
"""

params.i = 'input'
params.x = 'fasta'
params.db = "camitax/db"
params.dada2_db = 'rdp'


workflow {

    files = Channel
        .fromPath("${params.i}/*.${params.x}")
        .map { file -> [file.baseName, file] }
    db = Channel.fromPath( "${params.db}", type: 'dir').first()
    mash(db, files)
    prodigal(files)
    centrifuge(db, prodigal.out.fna_centrifuge)
    kaiju(db, prodigal.out.faa_kaiju)
    checkm(db, files)
    dada2_lineage = dada2(db, checkm.out.rRNA_fasta)
    
    id_collection = mash.out.mash_ANImax    .join(
    mash.out.mash_ANImax_taxIDs)            .join(
    dada2.out.dada2_lineage)                .join(
    centrifuge.out.centrifuge_taxIDs)       .join(
    kaiju.out.kaiju_taxIDs)                 .join(
    checkm.out.checkm_lineage)              .join(
    prodigal.out.genes_cnt)
    
    taxonomy(db, id_collection)
    summary(taxonomy.out.camitax_summaries.collect())
}

/********
 * (A) Genome-distance based assignment
 */
process mash {
    tag "${id}"

    input:
    file db
    tuple val(id), file(genome)

    output:
    tuple val(id), path("${id}.mash.ANImax.txt"), emit: mash_ANImax
    tuple val(id), path("${id}.mash.taxIDs.txt"), emit: mash_ANImax_taxIDs

    script:
    mash_index = "${db}/mash/RefSeq_20180510.msh"

    """
    mash dist ${mash_index} ${genome} | sort -gk3 > ${id}.mash.sorted
    head -n 1 ${id}.mash.sorted | awk '{print 1-\$3}' > ${id}.mash.ANImax.txt
    awk '\$3 <= 0.05' ${id}.mash.sorted | cut -f1 -d'_' > ${id}.mash.taxIDs.txt
    """
}

/********
 * (B) 16S rRNA gene-based assignment, and
 * (D) Phylogenetic placement
 */
process checkm {
    tag "${id}"

    input:
    file db
    tuple val(id), file(genome)

    output:
    tuple val(id), path("${id}.ssu.fna"), emit: rRNA_fasta
    tuple val(id), path("${id}.checkm.tsv"), emit: checkm_lineage

    script:
    checkm_db = "${db}/checkm/"

    """
    # Extract 16S rRNA gene sequences with Nhmmer
    checkm ssu_finder -t ${task.cpus} -x ${params.x} ${genome} . ssu_finder 2>&1
    ln -s ssu_finder/ssu.fna ${id}.ssu.fna

    # Phylogenetic placement onto reduced reference tree
    checkm lineage_wf -t ${task.cpus} --reduced_tree -x ${params.x} . checkm
    checkm qa -q -o 2 --tab_table checkm/lineage.ms checkm > ${id}.checkm.tsv
    """
}

process dada2 {
    tag "${id}"

    input:
    file db
    tuple val(id), file(genome)

    output:
    tuple val(id), path("${id}.dada2.txt"), emit: dada2_lineage

    script:
    if ( params.dada2_db == 'silva' ) {
        dada2_train_set = "${db}/dada2/silva_nr_v132_train_set.fa.gz"
        dada2_species_assignment = "${db}/dada2/silva_species_assignment_v132.fa.gz"
    } else if ( params.dada2_db == 'rdp' ) {
        dada2_train_set = "${db}/dada2/rdp_train_set_16.fa.gz"
        dada2_species_assignment = "${db}/dada2/rdp_species_assignment_16.fa.gz"
    } else
        error "Invalid database for dada2 specified: ${params.dada2_db}"

    """
    #!/usr/bin/env Rscript

    library(dada2)
    set.seed(42)

    seqs <- paste(Biostrings::readDNAStringSet("${genome}"))
    tt <- data.frame(Batman = "NA;NA;NA;NA;NA;NA;NA")
    try(tt <- addSpecies(assignTaxonomy(seqs, "${dada2_train_set}", minBoot=80), "${dada2_species_assignment}"))
    write.table(tt, "${id}.dada2.txt", quote=F, sep=";", row.names=F, col.names=F)
    """
}


/********
 * (C) Gene homology-based assignment(s)
 */
process prodigal {
    tag "${id}"
    publishDir "data/${id}"

    input:
    tuple val(id), file(genome)

    output:
    tuple val(id), path("${id}.genes.faa"), emit: faa_kaiju
    tuple val(id), path("${id}.genes.fna"), emit: fna_centrifuge
    tuple val(id), path("${id}.genes.cnt"), emit: genes_cnt

    """
    prodigal -a ${id}.genes.faa -d ${id}.genes.fna -i ${genome}
    grep -c "^>" ${id}.genes.faa > ${id}.genes.cnt
    """
}

process centrifuge {
    tag "${id}"

    input:
    file db
    tuple val(id), file(genes)

    output:
    tuple val(id), path("${id}.centrifuge.taxIDs.txt"), emit: centrifuge_taxIDs

    script:
    centrifuge_index = "${db}/centrifuge/proGenomes_20180510"

    """
    centrifuge -p ${task.cpus} -f -x ${centrifuge_index} ${genes} > ${id}.centrifuge.out
    awk '\$4+1 > 250' ${id}.centrifuge.out | cut -f3 > ${id}.centrifuge.taxIDs.txt
    """
}

process kaiju {
    tag "${id}"

    input:
    file db
    tuple val(id), file(genes)

    output:
    tuple val(id), path("${id}.kaiju.taxIDs.txt"), emit: kaiju_taxIDs

    script:
    kaiju_index = "${db}/kaiju/proGenomes_20180510.fmi"
    ncbi_nodes = "${db}/taxonomy/nodes_20180510.dmp"

    """
    kaiju -p -z ${task.cpus} -t ${ncbi_nodes} -f ${kaiju_index} -i ${genes} > ${id}.kaiju.out
    grep "^C" ${id}.kaiju.out | cut -f3 > ${id}.kaiju.taxIDs.txt
    """
}



process taxonomy {
    tag "${id}"

    input:
    file db
    tuple val(id), file(mash_ANImax), file(mash_taxIDs), file(dada2_lineage), file(centrifuge_taxIDs), file(kaiju_taxIDs), file(checkm_lineage), file(genes_cnt)

    output:
    tuple val(id), path("${id}.summary"), emit: camitax_summaries

    script:
    mash_ids = "${db}/mash/RefSeq_20180510.ids"
    ncbi_names = "${db}/taxonomy/names_20180510.dmp"
    ncbi_nodes = "${db}/taxonomy/nodes_20180510.dmp"

    """
    camitaxonomy.py --names ${ncbi_names} \
                    --nodes ${ncbi_nodes} \
                    --mash ${mash_taxIDs} \
                    --dada2 ${dada2_lineage} \
                    --centrifuge ${centrifuge_taxIDs} \
                    --kaiju ${kaiju_taxIDs} \
                    --checkm ${checkm_lineage} \
                    --known ${mash_ids} \
                    --animax ${mash_ANImax} \
                    --genes ${genes_cnt} \
                    ${id} > ${id}.summary
    """
}

process summary {
    tag 'The Final Countdown'
    publishDir 'data'

    input:
    file summaryList

    output:
    file "camitax.tsv"

    """
    for summary in ${summaryList}
    do
        head -n 1 \$summary > camitax.tsv
        break
    done
    for summary in ${summaryList}
    do
        tail -n +2 \$summary >> camitax.tsv
    done
    """
}

