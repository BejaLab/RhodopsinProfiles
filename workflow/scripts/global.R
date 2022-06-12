
library(dplyr)
library(tidyr)
library(seqinr)

with(snakemake@input, {
    usearch_file <<- usearch
    clstr_file   <<- clstr
    fasta_file   <<- fasta
    cdhit_file   <<- cdhit
})
with(snakemake@output, {
    similar_file   <<- similar
    different_file <<- different
})

read.cdhit.clstr <- function(fname) {
    read.table(fname, sep = "\t", comment.char = "", quote = "", fill = T, stringsAsFactors = F, col.names = c("Col1", "Col2")) %>%
        separate(Col1, into = c("Seq.Num", "Cluster"), sep = " ", fill = "right") %>%
        fill(Cluster) %>%
        filter(!grepl(">", Seq.Num)) %>%
        separate(Col2, into = c("Seq.Len", "Col2"), sep = "aa, >") %>%
        extract(Col2, into = c("Seq.Name", "Is.Representative", "Identity"), regex = "(.*?)[.]{3} ([*]|at) ?(.*)") %>%
        mutate(Is.Representative = Is.Representative == "*", Identity = ifelse(Is.Representative, "100%", Identity)) %>%
        group_by(Cluster) %>%
        mutate(Representative = Seq.Name[which(Is.Representative)]) %>%
        mutate(Identity = as.numeric(sub("%", "", Identity))) %>%
        ungroup
}

tsv <- read.table(usearch_file, col.names = c("Ref", "ID", "id"), sep = "\t", comment.char = "") %>%
    mutate(ID = sub(" .*", "", ID)) %>%
    mutate(Is.Ref = ID %in% Ref) %>%
    filter(Ref != ID) %>%
    mutate(similar = sprintf("%s [%s]", Ref, id))

clusters <- read.cdhit.clstr(clstr_file) %>%
    select(ID = Seq.Name, Cluster)

all_fasta <- read.fasta(fasta_file, as.string = T, seqtype = "AA") %>%
    {data.frame(ID = names(.), desc = sapply(., attr, "Annot"), full_sequence = unlist(.))} %>%
    extract(desc, into = c("Tax", "TaxID"), regex = "Tax=(.+?) TaxID=(\\d+)")

data <- read.fasta(cdhit_file, as.string = T, seqtype = "AA") %>%
    {data.frame(ID = names(.), domain_sequence = unlist(.))} %>%
    left_join(all_fasta, by = "ID") %>%
    left_join(tsv, by = "ID") %>%
    left_join(clusters, by = "ID") %>%
    group_by(ID, Cluster, Is.Ref, domain_sequence, full_sequence, Tax, TaxID) %>%
    summarize(similar = paste(na.omit(similar), collapse = "; ")) %>%
    ungroup %>%
    replace_na(list(similar = ""))
filter(data, similar != "") %>%
    write.table(similar_file, row.names = F, col.names = T, sep = "\t", quote = F)
filter(data, similar == "") %>%
    write.table(different_file, row.names = F, col.names = T, sep = "\t", quote = F)

