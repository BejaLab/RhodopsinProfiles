== Rhodopsin domain profiles ==

* To create the profile: `snakemake -c{cores} profile --config name={group}`
* To trim sequences in the database folder: `snakemake -c{cores} domain --config name={group}`

`{group}` corresponds to a subfolder in `input/` with two files expected: `references.fasta` and `outliers.txt`
