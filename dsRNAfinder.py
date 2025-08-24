import re

import pandas as pd
import argparse

def parse_gff(file_path):
    cds_records = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'CDS':
                continue

            seqid = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            cds_records.append((seqid, start, end, strand))

    return pd.DataFrame(cds_records, columns=['seqid', 'start', 'end', 'strand'])

def compute_distances(df):
    distances = []

    for seqid, group in df.groupby('seqid'):
        pos_strand = group[group['strand'] == '+']
        neg_strand = group[group['strand'] == '-']

        for _, plus_row in pos_strand.iterrows():
            for _, minus_row in neg_strand.iterrows():
                cds1_start, cds1_end = plus_row['start'], plus_row['end']
                cds2_start, cds2_end = minus_row['start'], minus_row['end']

                overlap = min(cds1_end, cds2_end) - max(cds1_start, cds2_start) + 1
                distances.append({
                    'seqid': seqid,
                    'cds1_start': cds1_start,
                    'cds1_end': cds1_end,
                    'cds1_strand': '+',
                    'cds2_start': cds2_start,
                    'cds2_end': cds2_end,
                    'cds2_strand': '-',
                    'distance': overlap
                })

    return pd.DataFrame(distances)

def main(input_gff, output_csv):
    gff = open(input_gff)
    for i in gff:
        if not re.match(r'^#', i):
            ls = i.rstrip().split("\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify CDS pairs on opposite strands and measure distance/overlap.")
    parser.add_argument("-i", "--input", required=True, help="Input GFF file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    main(args.input, args.output)
