#!/usr/bin/env nextflow

process cache_dada2_db {
    storeDir "${baseDir}/db2/dada2/"
	scratch = true

    output:
	file 'rdp_train_set_16.fa.gz' into rdp_train_set
    file 'rdp_species_assignment_16.fa.gz' into rdp_species_assignment

    """
	wget https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz
    wget https://zenodo.org/record/801828/files/rdp_species_assignment_16.fa.gz
    """
}
rdp_train_set.println { "Download complete: " + it }
rdp_species_assignment.println { "Download complete: " + it }


process cache_checkm_db {
    storeDir "${baseDir}/db2/checkm/"
	scratch = true

    output:
	file '*' into checkm_data

    """
	wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    tar xfv checkm_data_*.tar.gz && rm *.tar.gz
    """
}
checkm_data.println { "Download complete: " + it }
