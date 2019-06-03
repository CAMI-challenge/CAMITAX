#!/usr/bin/env nextflow

params.db = 'camitax/db'

process cache_taxonomy {
    storeDir "${params.db}"
    // scratch = true

    output:
    file 'taxonomy/names_20190108.dmp'
    file 'taxonomy/nodes_20190108.dmp'

    """
    wget https://zenodo.org/record/2597959/files/taxonomy_20190108.tar.gz
    tar xfv taxonomy_*.tar.gz && rm taxonomy_*.tar.gz
    """
}

process cache_mash_db {
    storeDir "${params.db}"
    // scratch = true

    output:
    file 'mash/RefSeq_20190108.ids'
    file 'mash/RefSeq_20190108.msh'

    """
    wget https://zenodo.org/record/2597959/files/mash_20190108.tar.gz
    tar xfv mash_*.tar.gz && rm mash_*.tar.gz
    """
}

process cache_dada2_db {
    storeDir "${params.db}/dada2"
    // scratch = true

    output:
    file 'rdp_train_set_16.fa.gz'
    file 'rdp_species_assignment_16.fa.gz'
    file 'silva_nr_v132_train_set.fa.gz'
    file 'silva_species_assignment_v132.fa.gz'

    """
    wget https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz
    wget https://zenodo.org/record/801828/files/rdp_species_assignment_16.fa.gz
    wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz
    wget https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz
    """
}

process cache_centrifuge_db {
    storeDir "${params.db}"
    // scratch = true

    output:
    file 'centrifuge/proGenomes_20190108.1.cf'
    file 'centrifuge/proGenomes_20190108.2.cf'
    file 'centrifuge/proGenomes_20190108.3.cf'
    file 'centrifuge/proGenomes_20190108.4.cf'

    """
    wget https://zenodo.org/record/2597959/files/centrifuge_20190108.tar.gz
    tar xfv centrifuge_*.tar.gz && rm centrifuge_*.tar.gz
    """
}

process cache_kaiju_db {
    storeDir "${params.db}"
    // scratch = true

    output:
    file 'kaiju/proGenomes_20190108.fmi'

    """
    wget https://zenodo.org/record/2597959/files/kaiju_20190108.tar.gz
    tar xfv kaiju_*.tar.gz && rm kaiju_*.tar.gz
    """
}

process cache_checkm_db {
    storeDir "${params.db}/checkm"
    // scratch = true

    output:
    file 'distributions/cd_dist.txt'
    file 'distributions/gc_dist.txt'
    file 'distributions/td_dist.txt'
    file 'genome_tree/genome_tree.derep.txt'
    file 'genome_tree/genome_tree.metadata.tsv'
    file 'genome_tree/genome_tree.taxonomy.tsv'
    file 'genome_tree/genome_tree_full.refpkg/CONTENTS.json'
    file 'genome_tree/genome_tree_full.refpkg/genome_tree.fasta'
    file 'genome_tree/genome_tree_full.refpkg/genome_tree.log'
    file 'genome_tree/genome_tree_full.refpkg/genome_tree.tre'
    file 'genome_tree/genome_tree_full.refpkg/phylo_modelEcOyPk.json'
    file 'genome_tree/genome_tree_reduced.refpkg/CONTENTS.json'
    file 'genome_tree/genome_tree_reduced.refpkg/genome_tree.fasta'
    file 'genome_tree/genome_tree_reduced.refpkg/genome_tree.log'
    file 'genome_tree/genome_tree_reduced.refpkg/genome_tree.tre'
    file 'genome_tree/genome_tree_reduced.refpkg/phylo_modelJqWx6_.json'
    file 'genome_tree/missing_duplicate_genes_50.tsv'
    file 'genome_tree/missing_duplicate_genes_97.tsv'
    file 'hmms/checkm.hmm'
    file 'hmms/checkm.hmm.ssi'
    file 'hmms/phylo.hmm'
    file 'hmms/phylo.hmm.ssi'
    file 'hmms_ssu/createHMMs.py'
    file 'hmms_ssu/SSU_archaea.hmm'
    file 'hmms_ssu/SSU_bacteria.hmm'
    file 'hmms_ssu/SSU_euk.hmm'
    file 'img/img_metadata.tsv'
    file 'pfam/Pfam-A.hmm.dat'
    file 'pfam/tigrfam2pfam.tsv'
    file 'selected_marker_sets.tsv'
    file 'taxon_marker_sets.tsv'
    file 'test_data/637000110.fna'

    """
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    tar xfv checkm_*.tar.gz && rm checkm_*.tar.gz
    """
}
