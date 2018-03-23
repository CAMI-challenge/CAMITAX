# CAMITAX
Taxon labels for microbial genomes

## Requirements

- [Nextflow](https://www.nextflow.io/)
- [barrnap](https://github.com/tseemann/barrnap)
- [bedtools](https://github.com/arq5x/bedtools2)
- [R](https://www.r-project.org/) with [dada2](https://www.bioconductor.org/packages/release/bioc/html/dada2.html) and [Biostrings](https://www.bioconductor.org/packages/release/bioc/html/Biostrings.html)
- [Python 3](https://www.python.org/) with [gzip](https://docs.python.org/3/library/gzip.html)

Most of these can be easily installed with [conda](https://conda.io/miniconda.html); I'll add more detailed instructions soon.

## Citation

Key methods (re)implemented in CAMITAX were described in our CAMI manuscript. Please cite:
* Sczyrba, Hofmann, Belmann, et al. (2017). **Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software.** *Nature Methods*, 14, 11:1063–1071. doi:[10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458)
