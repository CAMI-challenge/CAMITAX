#!/usr/bin/env nextflow

process cache_dada2_db {
    storeDir "${baseDir}/db2/dada2/"
	scratch = true

    output:
	file 'rdp_train_set_16.fa.gz' into dada2_db
    file 'rdp_species_assignment_16.fa.gz' into dada2_db

    """
	wget https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz
    wget https://zenodo.org/record/801828/files/rdp_species_assignment_16.fa.gz
    """
}
dada2_db.println { "Download complete: " + it }


process cache_checkm_db {
    storeDir "${baseDir}/db2/checkm/"
	scratch = true

    output:
	file '*' into checkm_db

    """
	wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    tar xfv checkm_data_*.tar.gz && rm *.tar.gz
    """
}
checkm_db.println { "Download complete: " + it }
