## Yunfei Dai
## 08/16/2021

"""
This script takes a list of protein (WP_) accessions (version) as input, and outputs their detailed information including length, source organism, and CDD domain predictions.

Usage:
python make_shortlist_info.py [infile] [CDD table] [outfile] 

Note: the infile should only contain protein accessions (WP_xxxxxx.x)
"""

import pandas as pd
from sys import argv
from Bio import Entrez
Entrez.email = "dai.yun@northeastern.edu"

try:
    shortlist = argv[1]
except IndexError:
    shortlist = "/scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/test_prot.txt"

try:
    CDD_domains = = argv[2]
except IndexError:
    print("Please provide CDD info!")
    
try:
    outfile = argv[3]
except IndexError:
    outfile = "/scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/test_output.txt"

# Read CDD table as pandas dataframe
df_CDD_domains = pd.read_csv(CDD_domains, sep='\t')
df_CDD_domains.rename(columns={"Short name": "domain_name"}, inplace=True)
# Since one protein can have multiple records in CDD domain, the following two lines ensures that indexs are unique
df_CDD_domains = df_CDD_domains.groupby(df_CDD_domains.Query)['domain_name'].apply(",".join).reset_index()
df_CDD_domains.set_index('Query', inplace=True)
print(df_CDD_domains.head())

# Get protein info from ncbi using efetch
with open(shortlist, 'r') as prot_list:
    for line in prot_list:
        prot = line.strip()
        handle = Entrez.efetch(db="protein", id=prot, retmode="xml")
        summary = Entrez.read(handle) 
        length = summary[0]['GBSeq_length']
        accession = summary[0]['GBSeq_primary-accession']
        acc_version = summary[0]['GBSeq_accession-version']
        seqid = summary[0]['GBSeq_other-seqids']
        organism = summary[0]['GBSeq_organism']
        taxonomy = summary[0]['GBSeq_taxonomy']

        '''
        # further parse information in "GBSeq_comment"
        for key in summary[0].keys():
            if key == "GBSeq_comment":
                [print(i) for i in summary[0][key].split(';')]
            elif key == "GBSeq_feature-table":
                for record in summary[0][key]:
                    [print(key2, " : ", record[key2]) for key2 in record.keys()]
            else:
                print(key, " : ", summary[0][key])
        '''

