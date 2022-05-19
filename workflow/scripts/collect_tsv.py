#!/usr/bin/env python3

import csv
import re

tsv_files = snakemake.input
out_file = str(snakemake.output)

best_records = {}
for tsv_file in tsv_files:
    with open(tsv_file) as fd:
        reader = csv.reader(fd, delimiter = "\t")
        header = next(reader)
        for row in reader:
            record = zip(headers, row)
            record_id = record['record_id']
            evalue = record['full_E_value']
            if record_id not in best_records or best_records[record_id]['full_E_value'] > full_E_value:
                best_records[record_id] = record
assert len(header) > 0, "No records found"
with open(out_file, 'w') as fd:
    tsv = csv.writer(fd, delimiter = "\t")
    tsv.writerow(header)
    for record in best_records.values():
        tsv.writerow(record.values())