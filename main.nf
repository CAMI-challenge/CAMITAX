#!/usr/bin/env nextflow
camitax_version = 'v0.2.0'
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
params.db = "$baseDir/db/"
params.dada2_db = 'rdp'

Channel
    .fromPath("${params.i}/*.${params.x}")
    .map { file -> [file.baseName, file] }
    .into {
        input_genomes;
        mash_genomes;
        barrnap_genomes;
        bedtools_genomes;
        prodigal_genomes;
        checkm_genomes;
        summary_genomes
    }

println("Input genomes:")
input_genomes.println { it[0] + " - " + it[1] }
db = Channel.fromPath( "${params.db}", type: 'dir').first()


/********
 * (A) Genome-distance based assignment
 */
process mash {
    tag "${id}"

    publishDir "data/${id}"
    cpus = 1
    memory = '2 GB'

    input:
    file db
    set val(id), file(genome) from mash_genomes

    output:
    set id, "${id}.mash.ANImax.txt", "${id}.mash.taxIDs.txt" into mash_ANImax_taxIDs

    script:
    mash_index = "${db}/mash/RefSeq.msh"

    """
    mash dist ${mash_index} ${genome} | sort -gk3 > ${genome.baseName}.mash.sorted
    head -n 1 ${genome.baseName}.mash.sorted | awk '{print 1-\$3}' > ${genome.baseName}.mash.ANImax.txt
    awk 'NR == 1 {t=\$3*1.1}; \$3 <= t && \$3 <= 0.05' ${genome.baseName}.mash.sorted | cut -f2 -d',' > ${genome.baseName}.mash.taxIDs.txt
    """
}


/********
 * (A) 16S rRNA gene-based assignment, and
 * (D) Phylogenetic placement
 */
process checkm {
    tag "${id}"

    publishDir "data/${id}"
    cpus = 4
    memory '16 GB'

    input:
    file db
    set val(id), file(genome) from checkm_genomes

    output:
    set id, "${id}.ssu.fna" into rRNA_fasta
    set id, "${id}.checkm.tsv" into checkm_lineage

    script:
    checkm_db = "${db}/checkm/"

    """
    # Set CheckM root data location
    echo ${checkm_db} | checkm data setRoot ${checkm_db}

    # Extract 16S rRNA gene sequences with Nhmmer
    checkm ssu_finder -t ${task.cpus} -x ${params.x} ${genome} . ssu_finder 2>&1
    ln -s ssu_finder/ssu.fna ${id}.ssu.fna

    # Phylogenetic placement onto reduced reference tree
    checkm lineage_wf -t ${task.cpus} --reduced_tree -x ${params.x} . checkm
    checkm qa -o 2 --tab_table checkm/lineage.ms checkm > ${id}.checkm.tsv
    """
}

// TODO Multithreading?
process dada2 {
    tag "${id}"

    publishDir "data/${id}"
    cpus = 1
    memory = '16 GB'

    input:
    file db
    set val(id), file(fasta) from rRNA_fasta

    output:
    set id, "${id}.dada2.txt" into dada2_lineage

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

    seqs <- paste(Biostrings::readDNAStringSet("${fasta}"))
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
    cpus = 1
    memory = '1 GB'

    input:
    set val(id), file(genome) from prodigal_genomes

    output:
    set id, "${id}.genes.faa" into faa_kaiju
    set id, "${id}.genes.fna" into fna_centrifuge

    """
    prodigal -a ${id}.genes.faa -d ${id}.genes.fna -i ${genome}
    """
}

process centrifuge {
    tag "${id}"

    publishDir "data/${id}"
    cpus = 4
    memory '16 GB'

    input:
    file db
    set val(id), file(genes) from fna_centrifuge

    output:
    set id, "${id}.centrifuge.taxIDs.txt" into centrifuge_ids

    script:
    centrifuge_index = "${db}/centrifuge/proGenomes"

    """
    centrifuge -p ${task.cpus} -f -x ${centrifuge_index} ${genes} > ${id}.centrifuge.out
    awk '\$4+1 > 250' ${id}.centrifuge.out | cut -f3 > ${id}.centrifuge.taxIDs.txt
    """
}

process kaiju {
    tag "${id}"

    publishDir "data/${id}"
    cpus = 4
    memory '16 GB'

    input:
    file db
    set val(id), file(genes) from faa_kaiju

    output:
    set id, "${id}.kaiju.taxIDs.txt" into kaiju_ids

    script:
    kaiju_index = "${db}/kaiju/proGenomes.fmi"
    ncbi_nodes = "${db}/taxonomy/nodes.dmp"

    """
    kaiju -p -z ${task.cpus} -t ${ncbi_nodes} -f ${kaiju_index} -i ${genes} > ${id}.kaiju.out
    grep "^C" ${id}.kaiju.out | cut -f3 > ${id}.kaiju.taxIDs.txt
    """
}


/********
 * (E) Classification algorithm
 */
mash_ANImax_taxIDs .join(
dada2_lineage      .join(
centrifuge_taxIDs  .join(
kaiju_taxIDs       .join(
checkm_lineage )))) .set{ id_collection }

process summarize {
    tag "${id}"

    publishDir "data/${id}"
    cpus = 1
    memory = '1 GB'
    container = 'python'

    input:
    file db
    set val(id), file(mash_ANImax), file(mash_taxIDs), file(dada2_lineage), file(centrifuge_taxIDs), file(kaiju_taxIDs), file(checkm_lineage) from id_collection

    output:
    file "${id}.mash.taxIDs.summary"
    file "${id}.dada2.taxIDs.txt"
    file "${id}.dada2.taxIDs.summary"
    file "${id}.centrifuge.taxIDs.summary"
    file "${id}.kaiju.taxIDs.summary"
    file "${id}.checkm.taxIDs.txt"
    file "${id}.checkm.taxIDs.summary"
    file "${id}.taxIDs.summary"

    script:
    ncbi_names = "${db}/taxonomy/names.dmp"
    ncbi_nodes = "${db}/taxonomy/nodes.dmp"

    """
    # Convert CheckM and Dada2 lineages to NCBI Taxonomy IDs
    tail -n +2 ${checkm_lineage} | cut -f2 | sed 's/.__//' | sed 's/ (UID.*//' > ${id}.checkm.txt
    to_ncbi_taxonomy.py ${id}.checkm.txt ${ncbi_names} > ${id}.checkm.taxIDs.txt
    to_ncbi_taxonomy.py ${dada2_lineage} ${ncbi_names} > ${id}.dada2.taxIDs.txt

    # Compute LCA, iuLCA and RTLpath for each software
    taxonomy_tools.py ${mash_taxIDs} ${ncbi_nodes} > ${id}.mash.taxIDs.summary
    taxonomy_tools.py ${id}.dada2.taxIDs.txt ${ncbi_nodes} > ${id}.dada2.taxIDs.summary
    taxonomy_tools.py ${centrifuge_taxIDs} ${ncbi_nodes} > ${id}.centrifuge.taxIDs.summary
    taxonomy_tools.py ${kaiju_taxIDs} ${ncbi_nodes} > ${id}.kaiju.taxIDs.summary
    taxonomy_tools.py ${id}.checkm.taxIDs.txt ${ncbi_nodes} > ${id}.checkm.taxIDs.summary

    # Extract LCA (for Mash and CheckM) and iuLCA (for the remaining)
    grep -w "^LCA" ${id}.mash.taxIDs.summary | cut -f2 >> ${id}.taxIDs.txt
    grep -w "^iuLCA" ${id}.dada2.taxIDs.summary | cut -f2 >> ${id}.taxIDs.txt
    grep -w "^iuLCA" ${id}.centrifuge.taxIDs.summary | cut -f2 >> ${id}.taxIDs.txt
    grep -w "^iuLCA" ${id}.kaiju.taxIDs.summary | cut -f2 >> ${id}.taxIDs.txt
    grep -w "^LCA" ${id}.checkm.taxIDs.summary | cut -f2 >> ${id}.taxIDs.txt

    # Generate global summary, also conatining the final RTLpath
    taxonomy_tools.py ${id}.taxIDs.txt ${ncbi_nodes} > ${id}.taxIDs.summary
    """
}
