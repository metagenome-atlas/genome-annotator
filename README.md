# Genome Annotator

[![Snakemake](https://img.shields.io/badge/snakemake-≥5-brightgreen.svg)](https://snakemake.bitbucket.io)

This workflow is designed to annotate the genes a set of genomes (nucleotide fasta files) with InterProScan and MetaCyc pathways and other genome properties.

## Authors

* Silas Kieser (@silask)

## Dependencies

For now, InterProScan v5 and [genomeProperties](https://genome-properties.readthedocs.io/en/latest/calculating.html#local-analysis-method) should be installed. In future this installation will be automated.

## Usage

    snakemake --config input="path/to/genome.fasta"

or for a set of genomes:

    snakemake --config input="path/to/genomes"

For more detailed configuration see the `config.yaml`

# Cluster execution

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster qsub --jobs 100

or

    snakemake --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.
