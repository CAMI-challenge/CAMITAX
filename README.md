# CAMITAX: Taxon labels for microbial genomes

![The CAMITAX taxonomic assignment workflow](workflow.png "The CAMITAX taxonomic assignment workflow")

**The CAMITAX taxonomic assignment workflow.**
CAMITAX assigns one NCBI Taxonomy ID (taxID) to each input genome *G* by combining phylogenetic placement with 16S rRNA gene-, gene homology-, and genome distance-based taxonomic assignments.

## Requirements

All you need is [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/). This is the recommended way to run CAMITAX.

#### Plan B

If Docker is no option, try [Sigularity](https://singularity.lbl.gov/) instead. Please first consult the [Nextflow documentation](https://www.nextflow.io/docs/latest/singularity.html) and then adjust the [CAMITAX configuration](nextflow.config) accordingly. In a nutshell, disable Docker and enable Singularity by replacing `docker.enabled = true` with `singularity.enabled = true`, and you should be all set.

#### Plan C

As a last resort, you may run CAMITAX without software containers. However, this is not recommended and you have to install [all software dependencies](requirements.txt) by yourself. We suggest [Bioconda](https://bioconda.github.io/) for this tedious task, the following (untested) code snippet should work:
```
while read requirement; do conda install --yes --channel "bioconda" $requirement; done < requirements.txt
```
If you need any help or further guidance: Please [get in touch](https://github.com/abremges/CAMITAX/issues)!

## User Guide

### Installation

```
nextflow pull abremges/CAMITAX
```

**TODO** Download script for CAMITAX databases.

### Input

CAMITAX expects all input genomes in (genomic/nucleotide multi-)FASTA format.
If your input genomes are in the folder `input/` with file extension `.fasta`, please run:
```
nextflow run abremges/CAMITAX --i input --x fasta
```

### Output

CAMITAX outputs a tab-seperated file containing the individual taxon assignments.

**TODO** Format and describe CAMITAX output.


## Citation

A manuscript describing the full scope of CAMITAX is currently in preparation.
In the meantime, please cite the GitHub repository and/or the CAMI manuscript:
* Sczyrba, Hofmann, Belmann, et al. (2017). **Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software.** *Nature Methods*, 14, 11:1063–1071. doi:[10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458)
