#!/usr/bin/env nextflow
camitax_version = 'v0.1.0'
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
params.dada2_db = 'rdp'

Channel
    .fromPath("${params.i}/*.${params.x}")
    .into {
        input_genomes;
        mash_genomes;
        barrnap_genomes;
        bedtools_genomes;
        prodigal_genomes;
        checkm_genomes;
        summary_genomes
    }

input_genomes.collect().println { "Input genomes: " + it }
db = Channel.fromPath( "$baseDir/db/", type: 'dir').first()


process mash {
    publishDir 'data', mode: 'copy'
    cpus = 1
    memory = '2 GB'

    input:
    file db
    file genome from mash_genomes

    output:
    file "${genome.baseName}.mash.ANImax.txt"
    file "${genome.baseName}.mash.taxIDs.txt" into mash_ids

    script:
    mash_index = "${db}/mash/RefSeq.msh"

    """
    mash dist ${mash_index} ${genome} | sort -gk3 > ${genome.baseName}.mash.sorted
    head -n 1 ${genome.baseName}.mash.sorted | awk '{print 1-\$3}' > ${genome.baseName}.mash.ANImax.txt
    awk 'NR == 1 {t=\$3*1.1}; \$3 <= t && \$3 <= 0.05' ${genome.baseName}.mash.sorted | cut -f2 -d',' > ${genome.baseName}.mash.taxIDs.txt
    """
}


process barrnap {
    publishDir 'data', mode: 'copy'
    cpus = 1
    memory = '1 GB'

    input:
    file genome from barrnap_genomes

    output:
    file "${genome.baseName}.barrnap.gff" into barrnap_gff

    """
    barrnap --threads 1 ${genome} > ${genome.baseName}.barrnap.gff
    """
}


process bedtools {
    publishDir 'data', mode: 'copy'
    cpus = 1
    memory = '1 GB'

    input:
    file genome from bedtools_genomes
    file "${genome.baseName}.barrnap.gff" from barrnap_gff

    output:
    file "${genome.baseName}.16S_rRNA.fasta" into rRNA_fasta

    """
    bedtools getfasta -fi ${genome} -bed <(grep "product=16S ribosomal RNA" ${genome.baseName}.barrnap.gff)  > ${genome.baseName}.16S_rRNA.fasta
    """
}


process dada2 {
    publishDir 'data', mode: 'copy'
    cpus = 1
    memory = '16 GB'

    input:
    file db
    file fasta from rRNA_fasta

    output:
    file "${fasta.baseName}.dada2.txt" into dada2_lineage

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
    write.table(tt, "${fasta.baseName}.dada2.txt", quote=F, sep=";", row.names=F, col.names=F)
    """
}


process to_ncbi_taxonomy {
    publishDir 'data', mode: 'copy'
    cpus = 1
    memory = '1 GB'
    container = 'python'

    input:
    file db
    file lineage from dada2_lineage

    output:
    file "${lineage.baseName}.taxIDs.txt" into dada2_ids

    script:
    ncbi_names = "${db}/taxonomy/names.dmp"

    """
    to_ncbi_taxonomy.py ${lineage} ${ncbi_names} > ${lineage.baseName}.taxIDs.txt
    """
}


process prodigal {
    publishDir 'data', mode: 'copy'
    cpus = 1
    memory = '1 GB'

    input:
    file genome from prodigal_genomes

    output:
    file "${genome.baseName}.prodigal.faa" into faa_kaiju
    file "${genome.baseName}.prodigal.fna" into fna_centrifuge

    """
    prodigal -a ${genome.baseName}.prodigal.faa -d ${genome.baseName}.prodigal.fna -f gff -i ${genome} -o ${genome.baseName}.prodigal.gff
    """
}


process centrifuge {
    publishDir 'data', mode: 'copy'
    cpus = 4
    memory '16 GB'

    input:
    file db
    file genes from fna_centrifuge

    output:
    file "${genes.baseName}.centrifuge.taxIDs.txt" into centrifuge_ids

    script:
    centrifuge_index = "${db}/centrifuge/proGenomes"

    """
    centrifuge -p ${task.cpus} -f -x ${centrifuge_index} ${genes} > ${genes.baseName}.centrifuge.out
    awk '\$4+1 > 250' ${genes.baseName}.centrifuge.out | cut -f3 > ${genes.baseName}.centrifuge.taxIDs.txt
    """
}


process kaiju {
    publishDir 'data', mode: 'copy'
    cpus = 4
    memory '16 GB'

    input:
    file db
    file genes from faa_kaiju

    output:
    file "${genes.baseName}.kaiju.taxIDs.txt" into kaiju_ids

    script:
    kaiju_index = "${db}/kaiju/proGenomes.fmi"
    ncbi_nodes = "${db}/taxonomy/nodes.dmp"

    """
    kaiju -p -z ${task.cpus} -t ${ncbi_nodes} -f ${kaiju_index} -i ${genes} > ${genes.baseName}.kaiju.out
    grep "^C" ${genes.baseName}.kaiju.out | cut -f3 > ${genes.baseName}.kaiju.taxIDs.txt
    """
}


process checkm {
    publishDir 'data', mode: 'copy'
    cpus = 4
    memory '16 GB'

    input:
    file db
    file genome from checkm_genomes

    output:
    file "${genome.baseName}.checkm.tsv"

    script:
    checkm_db = "${db}/checkm/"

    """
    echo ${checkm_db} | checkm data setRoot ${checkm_db}
    checkm lineage_wf -t ${task.cpus} --reduced_tree -x ${params.x} . checkm_out
    checkm qa -o 2 --tab_table checkm_out/lineage.ms checkm_out > ${genome.baseName}.checkm.tsv
    """
}


process summarize {
    publishDir 'data', mode: 'copy'
    container = 'python'

    input:
    file db
    file genome from summary_genomes
    file "${genome.baseName}.mash.taxIDs.txt" from mash_ids
    file "${genome.baseName}.16S_rRNA.dada2.taxIDs.txt" from dada2_ids
    file "${genome.baseName}.prodigal.centrifuge.taxIDs.txt" from centrifuge_ids
    file "${genome.baseName}.prodigal.kaiju.taxIDs.txt" from kaiju_ids

    output:
    file "${genome.baseName}.mash.taxIDs.summary"
    file "${genome.baseName}.dada2.taxIDs.summary"
    file "${genome.baseName}.centrifuge.taxIDs.summary"
    file "${genome.baseName}.kaiju.taxIDs.summary"
    file "${genome.baseName}.taxIDs.summary"

    script:
    ncbi_nodes = "${db}/taxonomy/nodes.dmp"

    """
    taxonomy_tools.py ${genome.baseName}.mash.taxIDs.txt ${ncbi_nodes} > ${genome.baseName}.mash.taxIDs.summary
    taxonomy_tools.py ${genome.baseName}.16S_rRNA.dada2.taxIDs.txt ${ncbi_nodes} > ${genome.baseName}.dada2.taxIDs.summary
    taxonomy_tools.py ${genome.baseName}.prodigal.centrifuge.taxIDs.txt ${ncbi_nodes} > ${genome.baseName}.centrifuge.taxIDs.summary
    taxonomy_tools.py ${genome.baseName}.prodigal.kaiju.taxIDs.txt ${ncbi_nodes} > ${genome.baseName}.kaiju.taxIDs.summary

    grep -w "^LCA" ${genome.baseName}.mash.taxIDs.summary | cut -f2 >> ${genome.baseName}.taxIDs.txt
    grep -w "^iuLCA" ${genome.baseName}.dada2.taxIDs.summary | cut -f2 >> ${genome.baseName}.taxIDs.txt
    grep -w "^iuLCA" ${genome.baseName}.centrifuge.taxIDs.summary | cut -f2 >> ${genome.baseName}.taxIDs.txt
    grep -w "^iuLCA" ${genome.baseName}.kaiju.taxIDs.summary | cut -f2 >> ${genome.baseName}.taxIDs.txt

    taxonomy_tools.py ${genome.baseName}.taxIDs.txt ${ncbi_nodes} > ${genome.baseName}.taxIDs.summary
    """
}
