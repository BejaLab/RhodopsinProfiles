Rhodopsin domain profiles
=========================

* To create the profiles: `snakemake -c{cores} --use-conda profile`. Training data (`references.fasta`, `config.yaml` and `outliers.txt`) should be located in `input/{name}/`.
* To extract domain sequences: `snakemake -c{cores} --use-conda domain`. The fasta files with sequences to trim should be put in `queries/` and have `.fasta` extension.
