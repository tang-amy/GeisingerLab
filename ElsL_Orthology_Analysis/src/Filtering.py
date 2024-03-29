## Yunfei Dai
## 2021/06/23


"""
This script takes a list of protein fasta sequences as input, and filters these protein sequences based on their length and whether they contain any signal sequence. 

Predictions for signal peptides are generated from prediction tools including SignalIP5.0 and Phobius.
http://www.cbs.dtu.dk/services/SignalP/
https://phobius.sbc.su.se/

Predictions for functional domains are generated from CDD domain prediction tool.
https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml

Usage example:

Filtering performed on 08/24/2021
only kept proteins between 105 aa - 200 aa

python /home/dai.yun/GeisingerLab/ElsL_Orthology_Analysis/src/Filtering.py \
-i /scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/GT9TM664013_blast.fasta \
-p /scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/prediction_results_blastp/GT9TM664013_blast-SIGNALP-GP--output_protein_type.txt \
-n /scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/prediction_results_blastp/GT9TM664013_blast-SIGNALP-GN--output_protein_type.txt \
-b /scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/prediction_results_blastp/GT9TM664013_phobius.txt \
-c /scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/prediction_results_blastp/GT9TM664013_CDDhits_standard.csv \
-t /scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/prediction_results_blastp/GT9TM664013_TMHMM_result.txt \
-x /scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/prediction_results_blastp/GT9TM664013_predisi_GN.txt \
-y /scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/prediction_results_blastp/GT9TM664013_predisi_GP.txt \
--min 105 --max 200 \
-o /scratch/dai.yun/2021July_ElsL_PhylogeneticAnalysis/refseq_GT7JRFP8013_blast/shortlist_new/GT9TM664013_shortlist_prot_acc.txt
"""
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from Bio import SeqIO, Entrez
import re
from optparse import OptionParser
from sys import argv

options = OptionParser()
options.add_option("-i", "--infile", dest="infile",
                   help="input is fasta file containing rpotein sequences")
options.add_option("-p", "--sigIP_pos", dest="pos",
                   help="Signal prediction results (Gram-positive)")
options.add_option("-n", "--sigIP_neg", dest="neg",
                   help="Signal prediction results (Gram-negative)")
options.add_option("-b", "--phobius", dest="phob",
                   help="phobius prediction results")
options.add_option("-c", "--cdd", dest="cdd",
                   help="CDD domain prediction results (.csv)")
options.add_option("-t", "--tmhmm", dest="tmhmm",
                   help="TMHMM prediction results")
options.add_option("-x", "--predisin", dest="predisin",
                   help="predisi prediction results (Gram-negative)")
options.add_option("-y", "--predisip", dest="predisip",
                   help="predisi prediction results (Gram-positive)")
options.add_option("--min", dest="min_length", type='int', default=0,
                   help="minimum sequence length to keep")
options.add_option("--max", dest="max_length", type='int', default=2000,
                   help="maximum sequence length to keep")
options.add_option("-H", "--histogram", dest="histogram", default="off",
                   help="option to plot histogram, default is off")
options.add_option("-o", "--outfile", dest="outfile",
                   help="specify the output file directory")

opts, args = options.parse_args()
# 3) Did Batch Entrez to get the FASTA sequences from the WP_ id's
fasta_sequence = opts.infile
# 4) Used the FASTA sequences as input into SignalP web server --
# did Gram-negative setting first, then Gram-positive setting;
# downloaded the results for each
SignalIP_pos = opts.pos
SignalIP_neg = opts.neg
# 5)  Used the FASTA sequences as input into Phobius web server (short output mode)
Phobius = opts.phob
# 6) Also used the FASTA sequences as input into CDD search (batch-CD search).
# [This way we have info that lets us exclude known PG-binding domains ("PG_binding_"; "LysM")]
CDD_domains = opts.cdd
TMHMM = opts.tmhmm
predisi_GN = opts.predisip
predisi_GP = opts.predisin
min_length = opts.min_length
max_length = opts.max_length
plot_switch = opts.histogram  # Default is "off" - the script by default won't generate a histogram of protein lengths, unless user uses "on".
outfile = opts.outfile

df_phobius = pd.read_csv(Phobius, sep=r"\s+", skiprows=0, index_col='SEQUENCE_ID')
df_signalIP_pos = pd.read_csv(SignalIP_pos, sep='\t', skiprows=1, index_col='# ID')
df_signalIP_pos.index.names = ["ID"]
df_signalIP_neg = pd.read_csv(SignalIP_neg, sep='\t', skiprows=1, index_col='# ID')
df_signalIP_neg.index.names = ["ID"]
with open(CDD_domains, 'r') as CDD_infile:
    if "#Batch CD-search tool" in CDD_infile.read():
        df_CDD_domains = pd.read_csv(CDD_domains, skiprows=7, sep='\t')
    else:
        df_CDD_domains = pd.read_csv(CDD_domains, sep='\t')
df_CDD_domains.rename(columns={"Short name": "domain_name"}, inplace=True)

ID_list = df_phobius.index.tolist()  # this is the list of 1000 hits
# make a filtered subset that excludes any protein that has:
# -a signal peptide, TAT signal, or lipoprotein signal from 4
# -an SP or TM from 5
# -a PG_binding domain from 6
subset_phobius = []  # subset that does not have an SP or TM from phobius (836 hits)
for hit in ID_list:
    if str(df_phobius.TM[hit]) == "0" and str(df_phobius.SP[hit]) == "0":
        subset_phobius.append(hit)

subset_signalIP = []  # subset of subset_phobius that also does not have SP/TAT/LIPO from SignalIP Prediction (823 hits)
for hit in subset_phobius:
    ID = hit.replace('|', '_')
    if (df_signalIP_pos.Prediction[ID] == "OTHER" and float(df_signalIP_pos.OTHER[ID]) > 0.75 and 
        df_signalIP_neg.Prediction[ID] == "OTHER" and float(df_signalIP_neg.OTHER[ID]) > 0.75):
        subset_signalIP.append(hit)

"""
def get_CDD_domains(query):
    #For a given sequence, return all predicted domains found in ElsL-refseq-select1000_CDDhitdata_standard.csv
    domains = []
    for entry in df_CDD_domains.index.tolist():
        if query in df_CDD_domains.Query[entry]:
            domains.append(df_CDD_domains.domain_name[entry])
    return domains
"""

df_CDD_domains['Query'] = df_CDD_domains['Query'].astype(str)
df_CDD_domains['domain_name'] = df_CDD_domains['domain_name'].astype(str)
df_CDD_domains = df_CDD_domains.groupby(df_CDD_domains.Query)['domain_name'].apply(",".join).reset_index()
df_CDD_domains.set_index('Query', inplace=True)

subset_CDD = []  # subset of subset_signalIP that does not contain PG_binding domain from CDD table (821 hits)
for hit in subset_signalIP:
    for entry in df_CDD_domains.index.tolist():
        if hit in entry:
            domain = df_CDD_domains.domain_name[entry]
            if 'PG_binding' not in domain:
                subset_CDD.append(hit)

# Additional filters
df_TMHMM = pd.read_csv(TMHMM, delimiter=r"\s+", skiprows=1,
                       names=['ID', 'len', 'ExpAA', 'First', 'PreHel', 'Topology'], index_col='ID', engine='python')
with open(predisi_GN, 'r') as fp_GN:
    search_word = "Truncation"
    if search_word in fp_GN.read():
        # regular expression to recognize either a tab or multiple spaces as delimiter
        df_predisi_GN = pd.read_csv(predisi_GN, delimiter=r"([ ]{2,})|(\t)", skiprows=7, index_col='FASTA-ID', engine='python')
    else:
        df_predisi_GN = pd.read_csv(predisi_GN, delimiter=r"([ ]{2,})|(\t)", skiprows=0, index_col='FASTA-ID', engine='python')

with open(predisi_GP, 'r') as fp_GP:
    search_word = "Truncation"
    if search_word in fp_GP.read():
        # regular expression to recognize either a tab or multiple spaces as delimiter
        df_predisi_GP = pd.read_csv(predisi_GP, delimiter=r"([ ]{2,})|(\t)", skiprows=7, index_col='FASTA-ID', engine='python')
    else:
        df_predisi_GP = pd.read_csv(predisi_GP, delimiter=r"([ ]{2,})|(\t)", skiprows=0, index_col='FASTA-ID', engine='python') 

subset_exclude_predisi_TMHMM = []
for hit in subset_CDD:
    if df_TMHMM.PreHel[hit] != "PredHel=1":
       for entry in df_predisi_GN.index.tolist():
           if hit in entry:
               if df_predisi_GN["Signal Peptide ?"][entry] == "N" and df_predisi_GP["Signal Peptide ?"][entry] == "N":
                   subset_exclude_predisi_TMHMM.append(hit)


# Read fasta file that contains protein sequences
fasta_sequences = SeqIO.parse(open(fasta_sequence), 'fasta')
seq_dict = {rec.id : rec.seq for rec in fasta_sequences}
all_seq_lengths = []

for seq_record in seq_dict:
    all_seq_lengths.append(len(seq_dict[seq_record]))

CDD_hit_length = []
for hit in subset_CDD:
    CDD_hit_length.append(len(seq_dict[hit]))

predisi_TMHMM_hit_length = []
for hit in subset_exclude_predisi_TMHMM:
    predisi_TMHMM_hit_length.append(len(seq_dict[hit]))

subset_size_excluded = []

"""
# for fasta input with description line like ">gi|490280925|ref|WP_004176841.1| L,D-transpeptidase [Nitrosospira lacus]"
for hit in subset_exclude_predisi_TMHMM:
    uid = hit.split("|")[1]
    protein_accession = hit.split("|")[3]
    size = len(seq_dict[hit])
    if size <= 225:
        subset_size_excluded.append(hit)
        for entry in df_CDD_domains.index.tolist():
            if hit in entry:
                species = re.findall(r"\[(.*?)\]", entry)[0]
                print(species)
                #outF.write(species)
                #outF.write("\n")
"""

# for fasta input with description line like ">WP_068936422.1 L,D-transpeptidase [Acinetobacter seifertii]"
outF = open(outfile, "a")
for hit in subset_exclude_predisi_TMHMM:
    protein_accession = hit
    size = len(seq_dict[hit])
    if min_length <= size <= max_length:
        subset_size_excluded.append(hit)
        outF.write(hit+"\n")
outF.close()

size_excluded_lengths = []
for hit in subset_size_excluded:
    length = len(seq_dict[hit])
    size_excluded_lengths.append(length)

if plot_switch == "on":
    plt.figure(figsize=(8,6))
    plt.hist(all_seq_lengths, bins='auto', alpha=0.5, label="all sequences")
    plt.hist(CDD_hit_length, bins='auto', alpha=0.5, label="subset_CDD_and_SigIP")
    plt.hist(predisi_TMHMM_hit_length, bins='auto', alpha=0.5, label="subset_CDD_SigIP_TMHMM_and_predisi")
    plt.hist(size_excluded_lengths, bins='auto', alpha=0.5, label="size_excluded (final shortlist)")

    plt.xlabel("Sequence Length")
    plt.ylabel("Count")
    plt.legend(loc="upper right")
    plt.show()
print("Total hits: " + str(len(all_seq_lengths)))
print("Hits passing all SP filters: " + str(len(subset_exclude_predisi_TMHMM)))
print("Hits shorter than max cutoff: " + str(len(size_excluded_lengths)))

