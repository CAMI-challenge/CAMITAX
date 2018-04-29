# CAMITAX: Taxon labels for microbial genomes

![The CAMITAX taxonomic assignment workflow](camitax.png "The CAMITAX taxonomic assignment workflow")

**The CAMITAX taxonomic assignment workflow.**
CAMITAX assigns one NCBI Taxonomy ID (taxID) to each input genome *G* by combining phylogenetic placement with 16S rRNA gene-, gene homology-, and genome distance-based taxonomic assignments.

## Requirements

All you need is [Nextflow](https://www.nextflow.io/) and [Docker](https://www.docker.com/). This is the recommended way to run CAMITAX.

#### Plan B

If Docker is no option, try [Sigularity](https://singularity.lbl.gov/) instead. Please consult the [Nextflow documentation](https://www.nextflow.io/docs/latest/singularity.html) and then adjust the [CAMITAX configuration](nextflow.config) accordingly. In a nutshell, disable Docker and enable Singularity by replacing `docker.enabled = true` with `singularity.enabled = true`.

#### Plan C

As a last resort, you may install all software dependencies manually and then run CAMITAX without software containers. We recommend [Bioconda](https://bioconda.github.io/) for this tedious task, good luck! Please [get in touch](https://github.com/abremges/CAMITAX/issues) if you need any help or further guidance.

## User Guide

### Installation

```
nextflow pull abremges/CAMITAX
```
That's all.

### Input

CAMITAX expects all input genomes in (genomic/nucleotide multi-)FASTA format.
If your input genomes are in the folder `input/` with file extension `.fasta`, please run:
```
nextflow run abremges/CAMITAX --i input --x fasta
```

### Output

CAMITAX outputs a tab-seperated file containing the individual taxon assignments. **TODO**


## Citation

A manuscript describing the full scope of CAMITAX is currently in preparation.
In the meantime, please cite the GitHub repository and/or the CAMI manuscript:
* Sczyrba, Hofmann, Belmann, et al. (2017). **Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software.** *Nature Methods*, 14, 11:1063–1071. doi:[10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458)
