Rhodopsin domain profiles
=========================

* To create the profile: `snakemake -c{cores} profile --configfile config/{group}.yaml`, where `{group}` is one of `7TM`, `8TM` or `HeR`. This assumes that the training set is located in `input/{group}/references.fasta` and a list of outliers is supplied in `input/{group}/outliers.txt`.
* To trim sequences: `snakemake -c{cores} domain --configfile config/{group}.yaml`. The sequences to trim are by default expected in `database/sequences.fasta`, otherwise specify with `--config sequences=/some/other/location/sequences.fasta`.
