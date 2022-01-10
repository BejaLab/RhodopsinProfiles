from Bio import SeqIO
import csv

a2m_file = str(snakemake.input["a2m"])
out_file = str(snakemake.output)

records = []
with open(a2m_file) as fh:
    a2m = SeqIO.parse(fh, 'fasta')
    ss_pred = next(a2m)
    ss_conf = next(a2m)
    cons = next(a2m)

    num_pos = len(cons)
    gap_nums = [0] * num_pos

    for record in a2m:
        records.append(record)
        for i in range(num_pos):
            gap_nums[i] += record.seq[i] == '-'

assert num_pos > 0, "Alignment is empty"

start = stop = None
for i in range(num_pos):
    gap_score = gap_nums[i] / num_records
    sec_struct = ss_pred.seq[i]
    if sec_struct == "H" and gap_score <= 0.1:
        start = i if start is None else start
        stop = i + 1

with open(out_file, "w") as out_fh:
    for record in records:
        record.seq = record.seq[start:stop]
        out_fh.write(record.format("fasta"))
