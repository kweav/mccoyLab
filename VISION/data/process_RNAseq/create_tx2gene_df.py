#!/usr/bin/env python3

'''
Usage: date;time ./create_tx2gene_df.py Mus_musculus.GRCm38.96.gtf
'''

import sys
import pandas as pd

'''
fields[2] needs to say transcript
then fields[8] can be mined for the gene id and the transcript id
gene_id "ENSMUSG00000102693"; gene_version "1"; transcript_id "ENSMUST00000193812"; transcript_version "1";
transcript_id needs to go to column1 and gene_id needs to go to column2
'''
transcript_ids = []
gene_ids = []
for line in open(sys.argv[1]):
    if line.startswith('#'):
        continue
    fields = line.strip('\r\n').split('\t')
    if fields[2] == 'transcript':
        subfields = fields[8].split(';')
        for field in subfields:
            split_it = field.split()
            if len(split_it) == 2:
                description, value = split_it[0], split_it[1]
                if description.strip() == "gene_id":
                    gene_id = value.replace('"','')
                elif description.strip() == "gene_version":
                    gene_version = value.replace('"','')
                elif description.strip() == "transcript_id":
                    transcript_id = value.replace('"','')
                elif description.strip() == "transcript_version":
                    transcript_version = value.replace('"','')
        if '{}.{}'.format(transcript_id, transcript_version) not in transcript_ids:
            transcript_ids.append('{}.{}'.format(transcript_id, transcript_version))
            gene_ids.append('{}.{}'.format(gene_id, gene_version))


transcript_to_gene_df = pd.DataFrame({'transcript_id': transcript_ids,
                                      'gene_id': gene_ids})

transcript_to_gene_df.to_csv('transcript_to_gene_df.Mus_musculus.GRCm38.96.csv', header=False, index=False)
