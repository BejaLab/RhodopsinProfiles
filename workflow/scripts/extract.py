#!/usr/bin/env python3

from Bio import SeqIO
from sys import stderr
import csv
import re

fas_file = str(snakemake.input["fasta"])
dom_file = str(snakemake.input["matches"])
tsv_file = str(snakemake.output)

offset          = snakemake.config["offset"]
extra_threshold = snakemake.config["extra_threshold"]
extra_offset    = snakemake.config["extra_offset"]
c_Evalue_threshold = snakemake.config["c_Evalue_threshold"]

domains = {}

# Parse the hmmsearch domain tabular output
with open(dom_file) as dom:
    for row in dom:
        if not row.startswith("#"):
            record_id, accession, tlen, query_name, accession, qlen,              \
                full_E_value, full_score, full_bias,                              \
                dom_num, dom_of, dom_c_Evalue, dom_i_Evalue, dom_score, dom_bias, \
                hmm_from, hmm_to, ali_from, ali_to, env_from, env_to,             \
                acc, description_of_target = re.split(" +", row, 23)
            if float(dom_c_Evalue) < c_Evalue_threshold:
                ali_from = int(ali_from)
                ali_to   = int(ali_to)
                qlen     = int(qlen)
                hmm_from = int(hmm_from)
                hmm_to   = int(hmm_to)
                hmm_left = qlen - hmm_to
                if record_id in domains and domains[record_id][1] < ali_to:
                    domains[record_id][1] = ali_to
                    domains[record_id][3] = hmm_left
                else:
                    domains[record_id] = [ ali_from, ali_to, hmm_from, hmm_left, full_E_value ]

# Parse the fasta file and extract the domains
with open(tsv_file, "w") as tsv_fh:
    tsv = csv.writer(tsv_fh, delimiter = "\t")
    tsv.writerow([ "record_id", "ali_from", "ali_left", "hmm_from", "hmm_left", "full_E_value", "domain_seq" ])
    with open(fas_file) as fas_fh:
        fas = SeqIO.parse(fas_fh, "fasta")
        for record in fas:
            if record.id in domains:
                ali_from, ali_to, hmm_from, hmm_left, full_E_value = domains[record.id]
                ali_left = len(record.seq) - ali_to

                offset_from = hmm_from + offset
                if hmm_from > extra_threshold:
                    offset_from += extra_offset

                offset_to = hmm_left + offset
                if hmm_left > extra_threshold:
                    offset_to += extra_offset

                trim_from = max(ali_from - offset_from, 0)
                trim_to   = ali_to + offset_to

                tsv.writerow([ record.id, str(ali_from), str(ali_left), str(hmm_from), str(hmm_left), full_E_value, str(record.seq[trim_from:trim_to]) ])
            else:
                tsv.writerow([ record.id, "", "", "", "", "", "" ])
                stderr.write("%s: domains not found\n" % record.id)
