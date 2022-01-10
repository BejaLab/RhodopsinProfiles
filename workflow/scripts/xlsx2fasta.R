library(readxl)
library(dplyr)
library(tidyr)

xlsx_file <- as.character(snakemake@input)
fasta_file <- as.character(snakemake@output)

sheets <- excel_sheets(xlsx_file)

lapply(sheets, function(sheet) read_xlsx(xlsx_file, sheet)) %>%
	setNames(sheets) %>%
	bind_rows(.id = "sheet") %>%
	filter(!is.na(ID)) %>%
	mutate(sprintf(">%s_%s\n%s", gsub(" ", "_", sheet), ID, Sequence)) %>%
	pull %>%
	cat(sep = "\n", file = fasta_file)
