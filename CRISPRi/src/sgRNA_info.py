## Yunfei Dai
## 04/15/2020


from optparse import OptionParser
import re
from bisect import bisect_left
import pandas as pd
from Bio import SeqIO
from src.sgRNA_finder import reverse_complement, locus_dic, locus_to_new

options = OptionParser()
options.add_option("-i", "--infile", dest="sgtable",
                   default='/Users/yunfei/GeisingerLab/CRISPRi/reference_files/CRISPRi_primer_info.tsv',
                   help="output file name)")
options.add_option("-t", "--tss", dest="tss",
                   default='/Users/yunfei/GeisingerLab/CRISPRi/reference_files/Kroger_TSS.txt',
                   help="input tss list (tab delimited)")
options.add_option("-r", "--reference", dest="reference",
                   default='/Users/yunfei/GeisingerLab/CRISPRi/reference_files/NZ_CP012004.fasta',
                   help="reference genome in .fasta format")
options.add_option("-g", "--genebank", dest="gbk",
                   default='/Users/yunfei/GeisingerLab/CRISPRi/reference_files/NZ_CP012004.gbk',
                   help='genome genebank file')
options.add_option("-o", "--outfile", dest="outfile",
                   default='/Users/yunfei/GeisingerLab/CRISPRi/results/sgRNA_info.tsv',
                   help="output file name)")


opts, args = options.parse_args()
sg_table = opts.sgtable
tss = opts.tss
reference = opts.reference
gbk = opts.gbk
outfile = opts.outfile
dic_locus = locus_dic(gbk)

df_sg = pd.read_csv(sg_table, sep='\t')
df_sg.drop(["Whole primer sequence ordered from IDT (5'-3')", "Distance of PAM from first base in ORF start codon",
            "Researcher name", "E. coli strain containing cloned construct", "Construct Notes",
            "phenotype by LB agar plating", "phenotype by LB broth growth", "other phenotype"], axis=1,
           inplace=True)
df_sg.rename(columns={"Primer Name (\"SGR-F\" primer containing the targeting sequence)": "Primer Name",
                      "gene is on which strand of genome (+ or -)": "Gene_Strand", "locus ID (ACX60_RSXXXX)": "locus ID",
                      "PAM site (on the + genome strand).  It is either NGG or CCN.  ": "PAM site",
                      "targeting sequence (5'-3')": "targeting sequence"}, inplace=True)
df_sg = df_sg[df_sg['locus ID'].notnull()]
df_sg.reset_index(inplace=True)
sg_index = df_sg.index.tolist()
sg_sequence = df_sg['targeting sequence'].tolist()
sg_Name = df_sg['Primer Name'].tolist()
sg_locus = df_sg['locus ID'].tolist()


df_TSS = pd.read_csv(tss, sep='\t', skipfooter=1, engine='python')
df_TSS.set_index('TSS coordinate', inplace=True)
tss_coordinate = df_TSS.index.tolist()
tss_plus = []
tss_minus = []
for t in tss_coordinate:
    if str(df_TSS['Strand'][t]) == '+':
        tss_plus.append(t)
    else:
        tss_minus.append(t)


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
                if gDNA.count(rc_seq) > 1:
                    print(seq, ' has more than one matches.')
                else:
                    type = 'CCN'
                    start_pos = gDNA.find(rc_seq)
                    end_pos = start_pos + len(rc_seq) - 1
                    N_pos = start_pos - 1
                    PAM = gDNA[N_pos - 2:N_pos + 1]
                    seed = seq[-12:]
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


def distance_finder(seq):
    global dic_locus
    info = pam_finder(reference, seq)
    i_plus = bisect_left(tss_plus, info[1])
    i_minus = bisect_left(tss_minus, info[2])
    match_plus = tss_plus[i_plus - 1]
    match_minus = tss_plus[i_minus]
    old_locus_plus = df_TSS['Locus_tag'][match_plus]
    new_locus_plus = locus_to_new(old_locus_plus, dic_locus)
    old_locus_minus = df_TSS['Locus_tag'][match_minus]
    new_locus_minus = locus_to_new(old_locus_minus, dic_locus)
    if abs(info[2] - match_plus) <= abs(match_minus - info[1]):
        if info[0] == 'NGG':
            return match_plus, '+', 'T', info[2] - match_plus, old_locus_plus, new_locus_plus, get_TSS_type(match_plus)
        else:
            return match_plus, '+', 'NT', match_plus - info[1], old_locus_plus, new_locus_plus, get_TSS_type(match_plus)
    else:
        if info[0] == 'NGG':
            return match_minus, '-', 'NT', info[2] - match_minus, old_locus_minus, new_locus_minus, get_TSS_type(match_minus)
        else:
            return match_minus, '-', 'T', match_minus - info[1], old_locus_minus, new_locus_minus, get_TSS_type(match_minus)


def get_TSS_type(pos):
    if df_TSS['Primary'][pos] == 1:
        tss_type = 'primary'
    elif df_TSS['Secondary'][pos] == 1:
        tss_type = 'secondary'
    elif df_TSS['Internal'][pos] == 1:
        tss_type = 'internal'
    elif df_TSS['Antisense'][pos] == 1:
        tss_type = 'antisense'
    elif df_TSS['Orphan'][pos] == 1:
        tss_type = 'orphan'
    else:
        tss_type = 'NA'
    return tss_type


def main():
    sgRNA_info = []
    for n in sg_index:
        name = sg_Name[n]
        seq = sg_sequence[n]
        locus_attempted = sg_locus[n]
        info_seq = pam_finder(reference, seq)
        start_pos = info_seq[1]
        end_pos = info_seq[2]
        pam = info_seq[3]
        seed = info_seq[4]
        Nearest_TSS = distance_finder(seq)[0]
        TSS_strand = distance_finder(seq)[1]
        target_strand = distance_finder(seq)[2]
        distance = distance_finder(seq)[3]
        old_tag = distance_finder(seq)[4]
        new_tag = distance_finder(seq)[5]
        type = distance_finder(seq)[6]
        sgRNA_info.append([name, seq, locus_attempted, start_pos, end_pos, pam, seed,
                           Nearest_TSS, TSS_strand, target_strand, distance, old_tag, new_tag, type])

    df_out = pd.DataFrame(sgRNA_info, columns=['Primer Name', 'Target Sequence (5\'-3\')', 'Attempted Target',
                                               'Start Position', 'End Position', 'PAM', 'Seed Number',
                                               'Nearest TSS Coordinate', 'TSS Strand', 'Target Strand',
                                               'Distance to TSS', 'Old Locus Tag', 'New Locus Tag', 'TSS Type'])
    df_out.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()