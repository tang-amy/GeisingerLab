import re
import pandas as pd

sg_table = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/CRISPRi_primer_info.tsv'
tss = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/Kroger_TSS.txt'
df_sg = pd.read_csv(sg_table, sep='\t')
df_sg.drop(["Whole primer sequence ordered from IDT (5'-3')", "Distance of PAM from first base in ORF start codon",
            "Researcher name", "E. coli strain containing cloned construct", "Construct Notes",
            "phenotype by LB agar plating", "phenotype by LB broth growth", "other phenotype"], axis=1, inplace=True)
df_sg.rename(columns={"gene is on which strand of genome (+ or -)": "Gene_Strand",
                      "PAM site (on the + genome strand).  It is either NGG or CCN.  ": "PAM site",
                      "targeting sequence (5'-3')": "targeting sequence"}, inplace=True)
df_sg = df_sg[df_sg['locus ID (ACX60_RSXXXX)'].notnull()]
df_tss = pd.read_csv(tss, sep='\t')
sg_sequence = df_sg['targeting sequence'].tolist()


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
                    start_pos = gDNA.find(reverse_complement(seq))
                    N_pos = start_pos - 1
                    PAM = gDNA[N_pos - 2:N_pos + 1]
                    seed = rc_seq[-12:]
                    rc_seed = reverse_complement(seed)
                    pattern_plus = ''.join([seed, '[ATGC]GG'])
                    pattern_minus = ''.join(['CC[ATGC]', rc_seed])
                    seed_number = len(re.findall(pattern_plus, gDNA)) + len(re.findall(pattern_minus, gDNA))
            else:
                print('Failed to locate short guide: ', seq)
                PAM = ''
                seed_number = 0
            return PAM, seed_number


def reverse_complement(seq):
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    seq = seq[::-1]
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)
