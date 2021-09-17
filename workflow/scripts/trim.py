from Bio import SeqIO
import csv

a2m_file = str(snakemake.input["a2m"])
sgc_file = str(snakemake.input["sgc"])
out_file = str(snakemake.output)

score = list()
with open(sgc_file) as fh:
    rd = csv.reader(fh, delimiter = "\t")
    for row in rd:
        if not row[0].startswith("+") and not row[0].startswith("|"):
            i = int(row[0])
            s = float(row[-1])
            score.insert(i, s)

num_pos = len(score)

assert num_pos > 0, "alignment is empty"

with open(out_file, "w") as out_fh:
    with open(a2m_file) as in_fh:
        a2m = SeqIO.parse(in_fh, "fasta")
        ss_pred = next(a2m)
        ss_conf = next(a2m)
        cons = next(a2m)

        assert len(ss_pred.seq) == num_pos, "ss_pred has a different length from trimal alignment scores: %d != %d" % (len(ss_pred.seq), num_pos)

        start = stop = None
        for i in range(num_pos):
            if ss_pred.seq[i] == "H" and score[i] > 0.9:
                start = i if start is None else start
                stop = i + 1

        for record in a2m:
            record.seq = record.seq[start:stop]
            out_fh.write(record.format("fasta"))
