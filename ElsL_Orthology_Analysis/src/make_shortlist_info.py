## Yunfei Dai
## 08/16/2021

"""
This script takes a list of protein (WP_) accessions (version) as input, and outputs their detailed information including length, source organism, and CDD domain predictions.

Usage:
python make_shortlist_info.py [infile] [CDD table] [outfile] 

Note: the infile should only contain protein accessions (WP_xxxxxx.x)
"""

import csv
import pandas as pd
from sys import argv
from Bio import Entrez
Entrez.email = "dai.yun@northeastern.edu"

try:
    shortlist = argv[1]
except IndexError:
    shortlist = "/scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/test_prot.txt"

try:
    CDD_domains = argv[2]
except IndexError:
    print("Please provide CDD info!")
    
try:
    outfile = argv[3]
except IndexError:
    outfile = "/scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/test_output.txt"

# Read CDD table as pandas dataframe
df_CDD_domains = pd.read_csv(CDD_domains, sep='\t', index_col=False)
df_CDD_domains.rename(columns={"Short name": "domain_name"}, inplace=True)


with open(shortlist, 'r') as prot_list:
    title= ['Accession', 'Accession_version', 'Length', 'SeqID', 'Organism', 'Superkingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
    record = []
    out_data = []
    for line in prot_list:
        # Get protein info from ncbi using efetch
        prot = line.strip()
        handle = Entrez.efetch(db="protein", id=prot, retmode="xml")
        summary = Entrez.read(handle) 
        length = summary[0]['GBSeq_length']
        accession = summary[0]['GBSeq_primary-accession']
        acc_version = summary[0]['GBSeq_accession-version']
        seqid = summary[0]['GBSeq_other-seqids']
        organism = summary[0]['GBSeq_organism']
        taxonomy = summary[0]['GBSeq_taxonomy'].split(";")
        superkingdom = taxonomy[0]
        phylum = taxonomy[1]
        tax_class = taxonomy[2]
        try:
            order = taxonomy[3]
        except IndexError:
            order = ""
        try:
            family = taxonomy[4]
        except IndexError:
            family = ""
        try:
            genus = taxonomy[5]
        except IndexError:
            genus = ""
        record = [accession, acc_version, length, ''.join(seqid), organism, superkingdom, phylum, tax_class, order, family, genus] 
        counter = 0
        # Get info from CDD prediction result
        for entry in df_CDD_domains.index.tolist():
            if prot in df_CDD_domains.Query[entry]:
                counter += 1
                CDD_title = ('CDD_record_'+str(counter))
                if CDD_title not in title:
                    title.append(CDD_title)
                evalue = df_CDD_domains.at[entry, 'E-Value']
                domain = df_CDD_domains.at[entry, 'domain_name']
                cdd_info = domain + '; eval='+str(evalue)
                record.append(cdd_info)
        out_data.append(record)
    # write data into output file
    df_output = pd.DataFrame(out_data)
    df_output.to_csv(outfile, header=title, sep='\t', index=False)
    
    
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

