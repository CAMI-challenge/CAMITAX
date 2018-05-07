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
        barrnap_genomes;
        bedtools_genomes;
        prodigal_genomes;
        mash_genomes
    }

input_genomes.collect().println { "Input genomes: " + it }
db = Channel.fromPath( "$baseDir/db/", type: 'dir').first()


process barrnap {
    publishDir 'data', mode: 'copy'

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
    maxForks 1

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
    container = 'python'

    input:
    file db
    file lineage from dada2_lineage

    output:
    file "${lineage.baseName}.ncbi.txt" into ncbi_lineage

    script:
    ncbi_names = "${db}/taxonomy/names.dmp"

    """
    to_ncbi_taxonomy.py ${lineage} ${ncbi_names} > ${lineage.baseName}.ncbi.txt
    """
}


// process prodigal {
//     publishDir 'data', mode: 'copy'
//
//     input:
//     file genome from prodigal_genomes
//
//     output:
//     file "${genome.baseName}.prodigal.faa" into prodigal_faa
//     file "${genome.baseName}.prodigal.fna" into prodigal_fna
//
//     """
//     prodigal -a ${genome.baseName}.prodigal.faa -d ${genome.baseName}.prodigal.fna -f gff -i ${genome} -o ${genome.baseName}.prodigal.gff
//     """
// }


// process checkm {
//     publishDir 'data', mode: 'copy'
//     maxForks 1
//
//     input:
//     file db
//     file genome from checkm_genomes
//     file "${genome.baseName}.prodigal.faa" from prodigal_faa
//
//     output:
//     file "${genome.baseName}.checkm.tsv"
//
//     script:
//     checkm_db = "${db}/checkm/"
//
//     """
//     echo ${checkm_db} | checkm data setRoot ${checkm_db}
//     checkm lineage_wf --reduced_tree --genes -x faa . checkm_out
//     checkm qa -o 2 --tab_table checkm_out/lineage.ms checkm_out > ${genome.baseName}.checkm.tsv
//     """
// }


// process centrifuge {
//     publishDir 'data', mode: 'copy'
//     maxForks 1
//
//     input:
//     file db
//     file genes from prodigal_fna
//
//     output:
//     file "${genes.baseName}.centrifuge.txt"
//
//     script:
//     centrifuge_index = "${db}/centrifuge/p_compressed"
//
//     """
//     centrifuge -f -x ${centrifuge_index} ${genes} > ${genes.baseName}.centrifuge.txt
//     """
// }


// process mash {
//     publishDir 'data', mode: 'copy'
//
//     input:
//     file db
//     file genome from mash_genomes
//
//     output:
//     file "${genome.baseName}.mash.ani"
//     file "${genome.baseName}.mash.tsv"
//
//     script:
//     mash_index = "${db}/mash/RefSeq.msh"
//
//     """
//     mash dist ${mash_index} ${genome} | sort -gk3 > ${genome.baseName}.mash.sorted
//     head -n 1 ${genome.baseName}.mash.sorted | awk '{print 1-\$3}' > ${genome.baseName}.mash.ani
//     awk 'NR == 1 {t=\$3*1.1}; \$3 <= t' ${genome.baseName}.mash.sorted > ${genome.baseName}.mash.tsv
//     """
// }
