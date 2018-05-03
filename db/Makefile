DIRGUARD = @mkdir -p $(@D)
DOWNLOAD = wget -O $@

all: \
dada2/rdp_species_assignment_16.fa.gz \
dada2/rdp_train_set_16.fa.gz \

dada2/rdp_species_assignment_16.fa.gz:
	$(DIRGUARD) && $(DOWNLOAD) https://zenodo.org/record/801828/files/rdp_species_assignment_16.fa.gz

dada2/rdp_train_set_16.fa.gz:
	$(DIRGUARD) && $(DOWNLOAD) https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz

#clean:
#	rm -rf dada2/
