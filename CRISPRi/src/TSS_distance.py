import re
import pandas as pd
from Bio import SeqIO
from src.sgRNA_finder import reverse_complement, locus_dic, locus_to_new

sg_table = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/CRISPRi_primer_info.tsv'
tss = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/Kroger_TSS.txt'
reference = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/NZ_CP012004.fasta'
gbk = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/NZ_CP012004.gbk'
dic_locus = locus_dic(gbk)

df_sg = pd.read_csv(sg_table, sep='\t')
df_sg.drop(["Whole primer sequence ordered from IDT (5'-3')", "Distance of PAM from first base in ORF start codon",
            "Researcher name", "E. coli strain containing cloned construct", "Construct Notes",
            "phenotype by LB agar plating", "phenotype by LB broth growth", "other phenotype"], axis=1,
           inplace=True)
df_sg.rename(columns={"gene is on which strand of genome (+ or -)": "Gene_Strand",
                      "PAM site (on the + genome strand).  It is either NGG or CCN.  ": "PAM site",
                      "targeting sequence (5'-3')": "targeting sequence"}, inplace=True)
df_sg = df_sg[df_sg['locus ID (ACX60_RSXXXX)'].notnull()]
sg_sequence = df_sg['targeting sequence'].tolist()


df_TSS = pd.read_csv(tss, sep='\t')
tss_index = df_TSS.index.tolist()[:-1]
tss_tag = [i for i in df_TSS["Locus_tag"].tolist()[:-1]]
tss_coordinate = [int(i) for i in df_TSS["TSS coordinate"].tolist()[:-1]]
tss_strand = df_TSS["Strand"][:-1].tolist()
tss_primary = [int(i) for i in df_TSS["Primary"].tolist()[:-1]]
tss_secondary = [int(i) for i in df_TSS["Secondary"].tolist()[:-1]]
tss_internal = [int(i) for i in df_TSS["Internal"].tolist()[:-1]]
tss_antisense = [int(i) for i in df_TSS["Antisense"].tolist()[:-1]]
tss_orphan = [int(i) for i in df_TSS["Orphan"].tolist()[:-1]]


def pam_finder(reference, seq):
    with open(reference, "r") as f:
        for record in SeqIO.parse(f, "fasta"):  # read sequence(s) from fasta
            gDNA = str(record.seq)
            rc_seq = reverse_complement(seq)
            if seq in gDNA:
                # Short guide is exact match, sequence preceding NGG (as on + strand)
                if gDNA.count(seq) > 1:
                    print(seq, ' has more than one matches.')
                else:
                    type = 'NGG'
                    start_pos = gDNA.find(seq)
                    end_pos = start_pos + len(seq) - 1
                    N_pos = start_pos + len(seq)
                    PAM = gDNA[N_pos:N_pos + 3]
                    seed = seq[-12:]
                    rc_seed = reverse_complement(seed)
                    pattern_plus = ''.join([seed, '[ATGC]GG'])
                    pattern_minus = ''.join(['CC[ATGC]', rc_seed])
                    seed_number = len(re.findall(pattern_plus, gDNA)) + len(re.findall(pattern_minus, gDNA))
            elif rc_seq in gDNA:
                # Short guide is reverse complementary match, sequence following CCN (as on + strand)
                if gDNA.count(seq) > 1:
                    print(seq, ' has more than one matches.')
                else:
                    type = 'CCN'
                    start_pos = gDNA.find(rc_seq)
                    end_pos = start_pos + len(rc_seq) - 1
                    N_pos = start_pos - 1
                    PAM = gDNA[N_pos - 2:N_pos + 1]
                    seed = rc_seq[-12:]
                    rc_seed = reverse_complement(seed)
                    pattern_plus = ''.join([seed, '[ATGC]GG'])
                    pattern_minus = ''.join(['CC[ATGC]', rc_seed])
                    seed_number = len(re.findall(pattern_plus, gDNA)) + len(re.findall(pattern_minus, gDNA))
            else:
                print('Failed to locate short guide: ', seq)
                type = None
                PAM = ''
                start_pos = 0
                end_pos = 0
                seed_number = 0
            return type, start_pos, end_pos, PAM, seed_number


def tss_finder(seq):
    info = pam_finder(reference, seq)
    if info[0] == 'NGG':
        pass
    elif info[0] == 'CCN':
        pass
    else:
        pass




