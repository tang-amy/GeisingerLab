## Yunfei Dai
## 03/27/2020

"""

This script identifies potential sgRNAs that target regions in proximity to transcription start sites (TSSs).

The default range queried is -50 to +100 of the TSSs (available from Kroger et al. 2017).

The script also checks for specificity of sgRNAs, only unique target sequences are kept.

Input file is .bed generated with PAM_finder.py

The output is a tab delimited txt file with 5 columns.
    Column 1: ACX_60 locus tag (old as in Kroger TSS list)
    Column 2: ACX_RS60 locus tag (new)
    Column 3: protein id
    Column 4: TSS strand
    Column 5: TSS coordinate
    Column 6: pam coordinate
    Column 7: pam sequence
    Column 8: SGR start position
    Column 9: SGR end position
    Column 10: target strand (template or non-template)
    Column 11: Distance (tss to the beginning / end of SGR sequence, whichever is longer)
    Column 11: SGR sequence
    Column 12: Seed number (total occurance of the 12 bp region preceeding PAM plus "NGG" in the entire genome)

[mandatory arguments]
    -i input file (.bed)
    -t TSS list (.csv)
    -r reference genome file (.fasta)
    -g genome annotation file (.gbk)
    -o output file

[optional arguments]
    -u acceptable upstream range to tss (default is 50)
    -d acceptable downstream range to tss (default is 100)
    -l length of desired sgRNA sequences (default is 20)

"""
from optparse import OptionParser
import os
import pandas as pd
from Bio import SeqIO
import timeit
import re

options = OptionParser()
options.add_option("-i", "--infile", dest="infile",
                   help="input list of PAM in .bed format)")
options.add_option("-t", "--tss", dest="tss",
                   help="input tss list (tab delimited)")
options.add_option("-r", "--reference", dest="reference",
                   help="reference genome in .fasta format")
options.add_option("-g", "--genebank", dest="gbk",
                   help='genome genebank file')
options.add_option("-u", "--up_range", dest="up_range", default="50")
options.add_option("-d", "--down_range", dest="down_range", default="100")
options.add_option("-l", "--sg_length", dest="sg_length", default="20")
options.add_option("-o", "--output_name", dest="outfile")
options.add_option("-n", "--TSS_list_type", dest="TSS_type", default="Kroger")

gbk = ''
TSS_type = ''

# Reads in pam_list.bed and TSS list.
# Identify pam sequences in proximity to at least one TSS.
# Output shortlist as "sg_candidate_[infile].bed"


def tss_pam(pam, tss, up_range, down_range, shortlist, sg_length):
    df_PAM = pd.read_csv(pam, sep='\t')
    df_TSS = pd.read_csv(tss, sep='\t')
    pam_index = df_PAM.index.tolist()  # index for all pam sequences
    all_pam = df_PAM["Start_pos"].tolist()  # start_pos of pam sequences
    pam_chrom = df_PAM["ID"].tolist()
    pam_end = df_PAM["End_pos"].tolist()
    pam_type = df_PAM["Type"].tolist()
    tss_index = df_TSS.index.tolist()[:-1]
    pam_coordinate = df_PAM["PAM_pos"].tolist()
    tss_coordinate = [int(i) for i in df_TSS["TSS coordinate"].tolist()[:-1]]
    tss_strand = df_TSS["Strand"][:-1].tolist()
    close_pam_index = []
    tss_plus = []
    tss_minus = []
    for t in tss_index:
        if tss_strand[t] == '+':
            tss_plus.append(tss_coordinate[t])
        else:
            tss_minus.append(tss_coordinate[t])
    for m in pam_index:
        start_to_search = all_pam[m]
        v = binary_search(tss_plus, 0, len(tss_plus) - 1, start_to_search, up_range, down_range, sg_length)
        w = binary_search(tss_minus, 0, len(tss_minus) - 1, start_to_search, down_range, up_range, sg_length)
        if v != -1:
            close_pam_index.append([m, '+'])
        if w != -1:
            close_pam_index.append([m, '-'])

    short_list = []
    for i in close_pam_index:
        chrom_id = pam_chrom[i[0]]
        pam_start = all_pam[i[0]]
        pam_end_pos = pam_end[i[0]]
        pam_t = pam_type[i[0]]
        coordinate = pam_coordinate[i[0]]
        strand_in_range = i[1]
        short_list.append([chrom_id, pam_t, pam_start, pam_end_pos, coordinate, strand_in_range])
    df = pd.DataFrame(short_list, columns=['ID', "Type", "Start_pos", "End_pos", "PAM_pos", "Strand_in_range"])
    print("Found " + str(len(short_list)) + " PAM sequences in proximity to at least one tss" + "\n")

    df.to_csv(shortlist, sep='\t', index=False, header=False)
    print("PAM shortlist saved as: " + shortlist + "\n")
    return shortlist  # returns path of the shortlist of candidates


def sgRNA_finder(sg_candidate, tss, reference, up_range, down_range, sg_length, sg_outfile):
    global gbk
    global TSS_type
    df_candidate = pd.read_csv(sg_candidate, sep='\t', names=['ID', "Type", "Start_pos", "End_pos", "PAM_pos", "Strand_in_range"])
    candidate_start = df_candidate["Start_pos"].tolist()  # start_pos of pam sequences
    candidate_end = df_candidate["End_pos"].tolist()
    candidate_type = df_candidate["Type"].tolist()
    candidate_coordinate = df_candidate["PAM_pos"].tolist()
    candidate_in_range_strand = df_candidate["Strand_in_range"].tolist()
    pam_index = df_candidate.index.tolist()
    t2 = timeit.default_timer()
    df_TSS = pd.read_csv(tss, sep='\t', skipfooter=1)
    tss_index = df_TSS.index.tolist()
    tss_coordinate = [int(i) for i in df_TSS["TSS coordinate"].tolist()]
    tss_strand = df_TSS["Strand"].tolist()
    old_locus = df_TSS["Locus tag"].tolist()
    dic_locus = locus_dic(gbk)
    if TSS_type == 'Kroger':
        dic_protein = pro_dic(gbk, 'old')
    elif TSS_type == 'Prados':
        dic_protein = pro_dic(gbk, 'new')
        TSS_cluster_ID = df_TSS['TSS cluster ID'].tolist()
        distance_ORF = df_TSS['Distance to ORF'].tolist()
        inside_ORF = df_TSS['is_inside_ORF'].tolist()

    # screen seed regions to avoid off-target effect (12 bp region next to NGG)
    # generates a dictionary {index: sequence}
    with open(reference, "r") as f:
        for record in SeqIO.parse(f, "fasta"):  # read sequence(s) from fasta
            gDNA = str(record.seq)
            candidate_seq = {}
            counter = 0
            for n in pam_index:
                sequence = str(record.seq[candidate_start[n]:candidate_end[n]])  # sgRNA sequence (default 20 bp)
                rc_sequence = reverse_complement(sequence)
                if candidate_type[n] == 'NGG':  # if NGG, keep original sequence
                    seed = sequence[-12:]
                    rc_seed = reverse_complement(seed)
                    pattern_plus = ''.join([seed, '[ATGC]GG'])
                    pattern_minus = ''.join(['CC[ATGC]', rc_seed])
                    seed_number = len(re.findall(pattern_plus, gDNA)) + len(re.findall(pattern_minus, gDNA))
                    if seed_number == 1:
                        counter += 1
                    p_seq = str(record.seq[(candidate_coordinate[n] - 1):(candidate_coordinate[n] + 2)])
                    candidate_seq.update({n: [sequence, p_seq, seed_number]})

                else:  # if CCN, keep reverse complementary
                    seed = rc_sequence[-12:]
                    rc_seed = reverse_complement(seed)
                    pattern_plus = ''.join([seed, '[ATGC]GG'])
                    pattern_minus = ''.join(['CC[ATGC]', rc_seed])
                    seed_number = len(re.findall(pattern_plus, gDNA)) + len(re.findall(pattern_minus, gDNA))
                    if seed_number == 1:
                        counter += 1
                    p_seq = str(record.seq[(candidate_coordinate[n] - 3):candidate_coordinate[n]])
                    candidate_seq.update({n: [rc_sequence, p_seq, seed_number]})

            t3 = timeit.default_timer()
            print("Finished screening seed regions in: ", format((t3 - t2), '.2f'), " seconds." + "\n")
            print("Found ", counter, " records with seed numbers = 1 .\n")

            if TSS_type == 'Kroger':
                df_sgRNA = pd.DataFrame(columns=['Locus tag (new)', 'Locus tag (old)', 'Protein ID', 'TSS strand',
                                                 'TSS-ORF distance', 'TSS coordinate', 'PAM coordinate', 'PAM sequence',
                                                 'SGR start', 'SGR end', 'Target strand', 'TSS-SGR distance',
                                                 'SGR sequence', 'Seed number', 'Additional TSS 1', 'Additional TSS 2',
                                                 'Additional TSS 3', 'Additional TSS 4', 'Additional TSS 5',
                                                 'Additional TSS 6'])
            elif TSS_type == 'Prados':
                df_sgRNA = pd.DataFrame(columns=['Locus tag', 'Protein ID', 'TSS strand', 'TSS coordinate',
                                                 'TSS-ORF distance', 'is_inside_ORF', 'PAM coordinate', 'PAM sequence',
                                                 'SGR start', 'SGR end', 'Target strand', 'TSS-SGR distance',
                                                 'SGR sequence', 'Seed number',  'TSS_cluster_ID', 'Additional TSS 1',
                                                 'Additional TSS 2', 'Additional TSS 3', 'Additional TSS 4'])
            index = 0
            for t in tss_index:
                for m in candidate_seq:
                    strand = tss_strand[t]
                    # locate sgRNA sequence to it's neighboring TSSs
                    if (-up_range < candidate_start[m] - tss_coordinate[
                        t] < down_range - sg_length and strand == '+') or \
                            (-down_range < candidate_start[m] - tss_coordinate[
                                t] < up_range - sg_length and strand == '-'):
                        if candidate_in_range_strand[m] == tss_strand[t]:
                            strand = tss_strand[t]
                            old_tag = old_locus[t]
                            if TSS_type == 'Kroger':
                                new_locus = locus_to_new(old_tag, dic_locus)
                            elif TSS_type == 'Prados':
                                cluster = TSS_cluster_ID[t]
                                inside = inside_ORF[t]
                                orf_dis = distance_ORF[t]
                            pro_id = locus_to_protein(old_tag, dic_protein)
                            tss_pos = tss_coordinate[t]
                            sg_start = candidate_start[m] + 1  # converted to 1-based index in output
                            sg_end = candidate_end[m]
                            pam_pos = candidate_coordinate[m]
                            p_seq = candidate_seq[m][1]
                            n_seed = candidate_seq[m][2]
                            if strand == '+':
                                t_strand = '+'
                                if pam_pos < sg_start < sg_end:
                                    p_type = "NT"
                                elif sg_start < sg_end < pam_pos:
                                    p_type = "T"
                            else:
                                t_strand = '-'
                                if pam_pos < sg_start < sg_end:
                                    p_type = "T"
                                elif sg_start < sg_end < pam_pos:
                                    p_type = "NT"
                            if sg_start <= tss_pos:
                                if strand == "+":
                                    distance = sg_start - tss_pos
                                else:
                                    distance = tss_pos - sg_start
                            else:
                                if strand == "+":
                                    distance = sg_end - tss_pos
                                else:
                                    distance = tss_pos - sg_end
                            sg_seq = ''.join(['5\'-', candidate_seq[m][0], '-3\''])
                            if TSS_type == 'Kroger':
                                df_sgRNA.loc[index] = [new_locus, old_tag, pro_id, t_strand, tss_pos, orf_dis,
                                                       pam_pos, p_seq, sg_start, sg_end, p_type, distance,
                                                       sg_seq, n_seed, cluster, '', '', '', '', '', '']
                            elif TSS_type == 'Prados':
                                if sg_seq not in df_sgRNA['SGR sequence'].tolist():
                                    df_sgRNA.loc[index] = [old_tag, pro_id, t_strand, tss_pos, orf_dis, inside,
                                                           pam_pos, p_seq, sg_start, sg_end, p_type, distance,
                                                           sg_seq, n_seed, cluster, '', '', '', '']
                                    index += 1
                                else:
                                    row_location = df_sgRNA.index[df_sgRNA['SGR sequence'] == sg_seq][0]
                                    old_distance = df_sgRNA['TSS-ORF distance'][row_location]
                                    if abs(old_distance) <= abs(distance):
                                        additional = ','.join(map(str,[old_tag, tss_pos, t_strand, orf_dis, distance]))
                                        if df_sgRNA['Additional TSS 1'][row_location] == '':
                                            df_sgRNA.at[row_location, 'Additional TSS 1'] = additional
                                        elif df_sgRNA['Additional TSS 2'][row_location] == '':
                                            df_sgRNA.at[row_location, 'Additional TSS 2'] = additional
                                        elif df_sgRNA['Additional TSS 3'][row_location] == '':
                                            df_sgRNA.at[row_location, 'Additional TSS 3'] = additional
                                        else:
                                            df_sgRNA.at[row_location, 'Additional TSS 4'] = additional
                                    else:
                                        additional = ', '.join(map(str, [df_sgRNA['Locus tag'][row_location],
                                                               df_sgRNA['TSS coordinate'][row_location],
                                                               df_sgRNA['TSS strand'][row_location],
                                                               df_sgRNA['TSS-ORF distance'][row_location],
                                                               df_sgRNA['TSS-SGR distance'][row_location]]))
                                        if df_sgRNA['Additional TSS 1'][row_location] == '':
                                            df_sgRNA.loc[row_location] = [old_tag, pro_id, t_strand, tss_pos, orf_dis,
                                                                          inside, pam_pos, p_seq, sg_start, sg_end,
                                                                          p_type, distance,sg_seq, n_seed, cluster,
                                                                          additional, '', '', '']
                                        elif df_sgRNA['Additional TSS 2'][row_location] == '':
                                            df_sgRNA.loc[row_location] = [old_tag, pro_id, t_strand, tss_pos, orf_dis,
                                                                          inside, pam_pos, p_seq, sg_start, sg_end,
                                                                          p_type, distance, sg_seq, n_seed, cluster,
                                                                          df_sgRNA['Additional TSS 1'][row_location],
                                                                          additional, '', '']
                                        elif df_sgRNA['Additional TSS 3'][row_location] == '':
                                            df_sgRNA.loc[row_location] = [old_tag, pro_id, t_strand, tss_pos, orf_dis,
                                                                          inside, pam_pos, p_seq, sg_start, sg_end,
                                                                          p_type, distance, sg_seq, n_seed, cluster,
                                                                          df_sgRNA['Additional TSS 1'][row_location],
                                                                          df_sgRNA['Additional TSS 2'][row_location],
                                                                          additional, '']
                                        else:
                                            df_sgRNA.loc[row_location] = [old_tag, pro_id, t_strand, tss_pos, orf_dis,
                                                                          inside, pam_pos, p_seq, sg_start, sg_end,
                                                                          p_type, distance, sg_seq, n_seed, cluster,
                                                                          df_sgRNA['Additional TSS 1'][row_location],
                                                                          df_sgRNA['Additional TSS 2'][row_location],
                                                                          df_sgRNA['Additional TSS 3'][row_location],
                                                                          additional]


# [Not executed] Assign ids to each unique sgRNA in sgRNA
# sg_start_list = []
# label = {}
# for entry in sgRNA:
#     sg_start_list.append(entry[8])
# temp = defaultdict(lambda: len(temp))
# for ele in sg_start_list:
#     label.update({ele: temp[ele]})
# for entry in sgRNA:
#     entry.insert(0, ''.join(['SGR_Ab', str(label[entry[8]]+1).zfill(4)]))

            df_sgRNA.to_csv(sg_outfile, sep='\t', index=False)
            t4 = timeit.default_timer()
            print("Finished recording sgRNA candidates: ", format((t4 - t3), '.2f'), " seconds.\n")
            print("Recorded ", len(candidate_seq), "sgRNA sequences.\n")
            print("sgRNA sequences saved as: ", sg_outfile, "\n")


def binary_search(arr_tss, l, r, p, up, down, length):
    if r >= l:
        mid = l + (r - l) // 2
        if -up <= p - arr_tss[mid] <= down - length:
            return mid
        elif p - arr_tss[mid] < -up:
            return binary_search(arr_tss, l, mid - 1, p, up, down, length)
        else:
            return binary_search(arr_tss, mid + 1, r, p, up, down, length)
    else:
        return -1


def reverse_complement(seq):
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    seq = seq[::-1]
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)


# Correspond old "ACX_60" locus tags to new "ACX_RS60" locus tags using .gbk genome file
# Ignore old tags that cannot be found in the gbk file
# Return a dictionary {old_tag: new_tag}
def locus_dic(annotation_gbk):
    recs = [rec for rec in SeqIO.parse(annotation_gbk, "genbank")]
    dic_locus = {}
    counter1 = 0
    counter2 = 0
    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == "gene"]
        for feat in feats:
            if 'old_locus_tag' in feat.qualifiers:
                old = ''.join((feat.qualifiers['old_locus_tag']))
                new = ''.join((feat.qualifiers['locus_tag']))
                counter1 += 1
            else:
                counter2 += 1
            dic_locus.update({old: new})

        dic_locus.update({old: new})
    print(counter1, "old tags converted to new tags. \n")
    print(counter2, "old tags failed to find match in the gbk reference, kept old tag. \n")
    return dic_locus


# Update a given old tag, return a new tag
def locus_to_new(locus_old, dic_locus):
    if dic_locus.get(locus_old) is None:
        locus_new = ''
    else:
        locus_new = dic_locus[locus_old]
    return locus_new


# Find protein_id based on old locus tag
# Return a dictionary {old_tag: protein_id}
def pro_dic(annotation_gbk, tag_version):
    recs = [rec for rec in SeqIO.parse(annotation_gbk, "genbank")]
    dic_protein = {}
    if tag_version == 'old':
        for rec in recs:
            feats = [feat for feat in rec.features if feat.type == "CDS"]
            for feat in feats:
                if 'old_locus_tag' in feat.qualifiers and 'protein_id' in feat.qualifiers:
                    locus = ''.join((feat.qualifiers['old_locus_tag']))
                    protein_id = ''.join((feat.qualifiers['protein_id']))
                dic_protein.update({locus: protein_id})
            return dic_protein
    elif tag_version == 'new':
        for rec in recs:
            feats = [feat for feat in rec.features if feat.type == "CDS"]
            for feat in feats:
                if 'locus_tag' in feat.qualifiers and 'protein_id' in feat.qualifiers:
                    locus = ''.join((feat.qualifiers['locus_tag']))
                    protein_id = ''.join((feat.qualifiers['protein_id']))
                dic_protein.update({locus: protein_id})
            return dic_protein


# Update a given old tag, return protein id
def locus_to_protein(locus, dic_protein):
    if dic_protein.get(locus) is None:
        protein_id = ''
    else:
        protein_id = dic_protein.get(locus)
    return protein_id


def main():
    opts, args = options.parse_args()
    pam = opts.infile
    pam_path = os.path.abspath(pam)
    pam_dir = os.path.dirname(pam_path)
    pam_filename = os.path.basename(pam_path)
    tss = opts.tss
    reference = opts.reference
    global gbk
    global TSS_type
    gbk = opts.gbk
    TSS_type = opts.TSS_type
    up_range = int(opts.up_range)
    down_range = int(opts.down_range)
    sg_length = int(opts.sg_length)
    sg_outfile = opts.outfile
    base_pam = os.path.splitext(pam_filename)[0]
    shortlist = pam_dir + '/' + 'PAM_shortlist_' + base_pam + '.bed'

    if not (os.path.isabs(sg_outfile)):
        sg_outfile = pam_dir + '/' + sg_outfile

    t0 = timeit.default_timer()
    shortlist = tss_pam(pam, tss, up_range, down_range, shortlist, sg_length)
    t1 = timeit.default_timer()
    print("finished PAM shortlist in: ", format(t1 - t0, '.2f'), " seconds" + "\n")
    sgRNA_finder(shortlist, tss, reference, up_range, down_range, sg_length, sg_outfile)


if __name__ == '__main__':
    main()
