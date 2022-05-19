#!/usr/bin/env python3

from Bio import SeqIO
from sys import stderr
from collections import defaultdict
import csv
import re

ref_hhr_file = str(snakemake.input['ref_hhr'])
a3m_file = str(snakemake.input['a3m'])
output_file = str(snakemake.output)
positions = snakemake.params['positions']

table = str.maketrans('', '', 'abcdefghijklmnopqrstuvwxyz')

hhm_positions = []

with open(ref_hhr_file) as hhr_fh:
	seq = defaultdict(str)
	pos = {}
	for line in hhr_fh:
		parts = line.rstrip().split()
		if len(parts) > 1 and parts[1] == "Consensus":
			qt, cons, start, cons_seq, end, total = parts
			seq[qt] += cons_seq
			pos[qt] = int(start) if qt not in pos else pos[qt]
	Q, T = seq['Q'], seq['T']
	i, j = pos['Q'], pos['T']
	for k in range(0, len(Q)):
		Q_no_gap = Q[k] != '-' and Q[k] != '.'
		T_no_gap = T[k] != '-' and T[k] != '.'
		if i in positions and Q_no_gap:
			assert T_no_gap, "Position %d aligns to gap" % i
			hhm_positions.append(j - 1)
		i += Q_no_gap
		j += T_no_gap
with open(a3m_file) as a3m_fh:
	with open(output_file, 'w') as out_fh:
		for record in SeqIO.parse(a3m_fh, "fasta"):
			matches = str(record.seq).translate(table)
			residues = []
			for pos in hhm_positions:
				residues.append(matches[pos])
			out_fh.write("%s\t%s\t%s\n" % (record.id, ','.join(residues), record.seq))
