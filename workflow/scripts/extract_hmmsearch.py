#!/usr/bin/env python3

from Bio import SeqIO
from sys import stderr
import csv
import re

fas_file = str(snakemake.input['fasta'])
dom_file = str(snakemake.input['matches'])
a2m_file = str(snakemake.input['a2m'])
tsv_file = str(snakemake.output)

config = snakemake.params['config']
name            = config['name']
offset_n        = config['offset_n']
offset_c        = config['offset_c']
extra_threshold = config['extra_threshold']
extra_offset    = config['extra_offset']
overlap         = config['overlap']
c_Evalue_threshold = config['c_Evalue_threshold']
hits_only = config['hits_only'] if 'hits_only' in config else False

positions = snakemake.params['positions']

domains = {}

hmm_len = None

def parse_hmm_len(line):
    """Parse the profile length line."""
    match = re.findall('\[M=(\d+)\]', line)
    return(int(match[-1]))

def parse_full_scores(fh):
    """Parse the full scores block.

    The block ends with an empty line.
    """
    scores = {}

    field_names = 'full_E_value', 'full_score', 'full_bias', 'best_E_value', 'best_score', 'best_bias', 'exp', 'N', 'Sequence'
    while (line := next(fh).strip()) != "":
        if not line.startswith('--') and not line.startswith('E-value'):
            fields = re.split(' +', line)
            assert len(fields) >= len(field_names), "Unexpected number of fields on line %s" % line
            vals = dict(zip(field_names, fields[0:len(field_names)]))
            scores[vals['Sequence']] = vals
    return(scores)

def parse_sequence_data(fh):
    """Parse a block for one input sequence.
    
    The block starts with ">>" and
    ends with last domain's last alignment block.
    """
    domains = []
    field_names = 'num', 'passed', 'score', 'bias', 'c_Evalue', 'i_Evalue', 'hmm_from', 'hmm_to', 'hmm_compl', 'ali_from', 'ali_to', 'ali_compl', 'env_from', 'env_to', 'env_compl', 'acc'

    no_domains = False

    # The scores sections starts immediately and ends with an empty line
    while (line := next(fh).strip()) != '' and not no_domains:
        no_domains = line.startswith("[No individual domains")
        if not line.startswith('#') and not line.startswith('--') and not no_domains:
            fields = re.split(' +', line)
            assert len(fields) == len(field_names), "Unexpected number of fields on line %s" % line
            vals = dict(zip(field_names, fields))
            vals['hmm_from'] = int(vals['hmm_from'])
            vals['hmm_to']   = int(vals['hmm_to'])
            vals['ali_from'] = int(vals['ali_from'])
            vals['ali_to']   = int(vals['ali_to'])
            vals['env_from'] = int(vals['env_from'])
            vals['env_to']   = int(vals['env_to'])
            vals['c_Evalue'] = float(vals['c_Evalue'])
            domains.append(vals)

    alignments_for_each_domain = next(fh) # skip this line

    i = -1
    while i < len(domains) - 1:
        line = next(fh).strip()
        if line.startswith('=='):
            i += 1
            domains[i]['RF'] = ''
            domains[i]['PP'] = ''
            domains[i]['hmm_seq'] = ''
            domains[i]['ali_seq'] = ''
            domains[i]['matches'] = ''
            domains[i]['hmm_right'] = -1
            while domains[i]["hmm_right"] < domains[i]['hmm_to']:
                block = parse_aln_block(fh)
                domains[i]['RF'] += block['RF']
                domains[i]['PP'] += block['PP']
                domains[i]['hmm_seq'] += block['hmm_seq']
                domains[i]['ali_seq'] += block['ali_seq']
                domains[i]['matches'] += block['matches']
                if 'hmm_left' not in domains[i]:
                    domains[i]['hmm_left'] = block['hmm_left']
                if 'ali_left' not in domains[i]:
                    domains[i]['ali_left'] = block['ali_left']
                domains[i]['hmm_right'] = block['hmm_right']
                domains[i]['ali_right'] = block['ali_right']
    return(domains)

def parse_aln_block(fh):
    """Parse profile-query alignment (sub)block.

    The block ends with an empty line.
    """
    vals = {}

    line = next(fh)
    fields = re.split(' +', line.strip(), 2)
    vals['RF'] = fields[0]

    line = next(fh)
    fields = re.split(" +", line.strip())
    vals['hmm_left']  = int(fields[-3])
    vals['hmm_seq']   = fields[-2]
    vals['hmm_right'] = int(fields[-1])

    line = next(fh)
    vals['matches'] = line.rstrip('\n\r')[-len(vals['hmm_seq']):]

    line = next(fh)
    fields = re.split(" +", line.strip())
    # Sequence line can have '-' for positions when alignment is all gaps (no sequence aligned)
    # Gap-only sequences are handled correctly downstream (gaps are skipped in position tracking)
    # Use ali_left/ali_right = None to indicate missing/invalid values
    ali_left_str = fields[-3]
    ali_right_str = fields[-1]
    vals['ali_left']  = int(ali_left_str) if ali_left_str != '-' else None
    vals['ali_seq']   = fields[-2]
    vals['ali_right'] = int(ali_right_str) if ali_right_str != '-' else None
    
    line = next(fh)
    fields = re.split(' +', line.strip(), 2)
    vals['PP'] = fields[0]

    line = next(fh).strip()
    assert line == '', "An empty line expected, got %s" % line

    return(vals)

# parse the hmm file to extract the Lys position
hmm_positions = []
with open(a2m_file) as fh:
    a2m = SeqIO.parse(fh, 'fasta')
    record = next(a2m)
    br_pos = 0
    hmm_pos = 0
    for res in record.seq:
        is_gap = res == '-'
        br_pos += not is_gap
        hmm_pos += not res.islower()
        if not is_gap and br_pos in positions:
            hmm_positions.append(hmm_pos)

# Parse the hmmsearch output
data = {}
with open(dom_file) as dom:
    full_scores = {}
    for line in dom:
        if line.startswith('Query:'):
            hmm_len = parse_hmm_len(line)
        elif line.startswith('Scores'):
            full_scores = parse_full_scores(dom)
        elif line.startswith('>>'):
            sequence_name = re.match('>>\s+(.+?)\s+', line)[1]
            data[sequence_name] = full_scores[sequence_name]
            data[sequence_name]['domains'] = parse_sequence_data(dom)

# Parse the fasta file and extract the domains
with open(fas_file) as fas_fh:
    with open(tsv_file, 'w') as tsv_fh:
        tsv = csv.writer(tsv_fh, delimiter = "\t")
        tsv.writerow([ 'record_id', 'profile', 'ali_from', 'ali_left', 'hmm_from', 'hmm_left', 'full_E_value', 'full_score', 'positions', 'positions_trim', 'positions_res', 'domain_seq' ])
        fas = SeqIO.parse(fas_fh, 'fasta')
        for record in fas:
            if record.id in data:
                blocks = []
                for domain in data[record.id]['domains']:
                    if domain['c_Evalue'] <= c_Evalue_threshold:
                        blocks.append(( domain['hmm_from'], domain['hmm_to'], domain['ali_from'], domain['ali_to'], domain['hmm_seq'], domain['ali_seq'] ))

                blocks.sort(key = lambda x: x[0])
                groups = []
                for i in range(len(blocks)):
                    hmm1_from, hmm1_to, ali1_from, ali1_to, *rest = blocks[i]
                    connected = False
                    for group in groups:
                        if i in group:
                            connected = True
                    if not connected:
                        for j in range(i + 1, len(blocks)):
                            hmm2_from, hmm2_to, ali2_from, ali2_to, *rest = blocks[j]
                            if hmm1_to - hmm2_from < overlap and ali1_to - ali2_from < overlap:
                                chained = False
                                for group in groups:
                                    if i in group or j in group:
                                        group.add(i)
                                        group.add(j)
                                        chained = True
                                if not chained:
                                    groups.append({i, j})
                                connected = True
                        if not connected:
                            groups.append({i})

                ranges = []
                for group in groups:
                    hmm_from = min([ blocks[i][0] for i in group ])
                    hmm_to   = max([ blocks[i][1] for i in group ])
                    ali_from = min([ blocks[i][2] for i in group ])
                    ali_to   = max([ blocks[i][3] for i in group ])
                    hmm_seqs = [ blocks[i][4] for i in group ]
                    ali_seqs = [ blocks[i][5] for i in group ]
                    ranges.append((hmm_from, hmm_to, ali_from, ali_to, hmm_seqs, ali_seqs ))

                ranges.sort(key = lambda x: x[2] - x[3])
                ranges_done = []
                for hmm_from, hmm_to, ali_from, ali_to, hmm_seqs, ali_seqs in ranges:
                    nested = False
                    for hmm0_from, hmm0_to, ali0_from, ali0_to in ranges_done:
                        if hmm0_from <= hmm_from <= hmm0_to and hmm0_from <= hmm_to <= hmm0_to:
                            nested = True
                    if not nested:
                        ali_positions      =   [0] * len(positions)
                        ali_positions_trim =   [0] * len(positions)
                        ali_positions_res  = ['-'] * len(positions)
                        for i in range(len(hmm_seqs)):
                            ali_pos = ali_from
                            hmm_pos = hmm_from
                            for j in range(len(hmm_seqs[i])):
                                hmm_res = hmm_seqs[i][j]
                                ali_res = ali_seqs[i][j]
                                if hmm_pos in hmm_positions and hmm_res != '.':
                                    pos_index = hmm_positions.index(hmm_pos)
                                    ali_positions[pos_index] = ali_pos
                                    ali_positions_res[pos_index] = ali_res
                                hmm_pos += hmm_res != '.'
                                ali_pos += ali_res != '-'

                        ali_left = len(record.seq) - ali_to
                        hmm_left = hmm_len - hmm_to

                        offset_from = offset_n

                        if hmm_from > extra_threshold:
                            offset_from += extra_offset

                        offset_to = offset_c
                        if hmm_left > extra_threshold:
                            offset_to += extra_offset

                        trim_from = max(ali_from - offset_from - hmm_from, 0)
                        for i in range(len(positions)):
                            ali_positions_trim[i] = ali_positions[i] - trim_from if ali_positions[i] > 0 else 0

                        trim_to   = ali_to + offset_to + hmm_left

                        new_seq = record.seq[trim_from:ali_from - 1].lower() + record.seq[ali_from - 1:ali_to] + record.seq[ali_to:trim_to].lower()
                        ali_positions_str      = ','.join(map(str, ali_positions))
                        ali_positions_trim_str = ','.join(map(str, ali_positions_trim))
                        ali_positions_res_str  = ','.join(ali_positions_res)
                        tsv.writerow([ record.id, name, str(ali_from), str(ali_left), str(hmm_from), str(hmm_left), data[record.id]['full_E_value'], data[record.id]['full_score'], ali_positions_str, ali_positions_trim_str, ali_positions_res_str, new_seq ])

            elif not hits_only:
                tsv.writerow([ record.id, '', '', '', '', '', '', '', '', '', '', '' ])
                stderr.write("%s: domains not found\n" % record.id)
