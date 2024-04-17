import csv
import re
from collections import defaultdict

tsv_files = snakemake.input
out_file = str(snakemake.output)

records = defaultdict(lambda: defaultdict(list))
for tsv_file in tsv_files:
    profile_records = defaultdict(list)
    with open(tsv_file) as fd:
        reader = csv.reader(fd, delimiter = "\t")
        header = next(reader)
        for row in reader:
            record = dict(zip(header, row))
            if record['profile'] != '':
                record_id = record['record_id']
                records[record_id][record['profile']].append(record)

assert len(header) > 0, "No records found"
with open(out_file, 'w') as fd:
    tsv = csv.writer(fd, delimiter = "\t")
    tsv.writerow(header)
    for record_id, record_data in records.items():
        scores = { profile: float(records[0]['full_score']) for profile, records in record_data.items() }
        best_profile = max(scores, key = scores.get)
        for record in record_data[best_profile]:
            tsv.writerow(record.values())
