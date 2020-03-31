## Yunfei Dai
## 03/27/2020

"""

This script identifies potential sgRNAs that target regions in proximity to transcription start sites (TSSs).

The default range queried is -50 to +100 of the TSSs (available from Kroger et al. 2017).

The script also checks for specificity of sgRNAs, only unique target sequences are kept.

Input file is .bed generated with PAM_finder.py

The output is a tab delimited txt file with 5 columns.
    Column 1: chromosome ID
    Column 2: "+" or "-" strand
    Column 3: sgRNA start position
    Column 4: sgRNA end position
    Column 5: sgRNA sequence
    Column 6: ACX60_ locus_tag (could be multiple for one sgRNA)

[mandatory options]
    -i input file (.bed)
    -t TSS list (.csv)
    -r reference genome file (.fasta)
    -g genome annotation file (.gbk)
    -o output file

[optional]
    -u acceptable upstream range to tss (default is 50)
    -d acceptable downstream range to tss (default is 100)

"""
from optparse import OptionParser
import os
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import timeit

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
options.add_option("-d", "--down_range", dest="down_range", default="100"),
options.add_option("-o", "--output_name", dest="outfile")

gkb = ""


# Reads in pam_list.bed and TSS list.
# Identify pam sequences in proximity to at least one TSS.
# Output shortlist as "sg_candidate_[infile].bed"

def tss_pam(pam, tss, up_range, down_range, shortlist):
    df_PAM = pd.read_csv(pam, sep='\t', names=["ID", "Strand", "Start_pos", "End_pos", "PAM_pos"])
    df_TSS = pd.read_csv(tss, sep='\t')
    pam_index = df_PAM.index.tolist()  # index for all pam sequences
    all_pam = df_PAM["Start_pos"].tolist()  # start_pos of pam sequences
    pam_chrom = df_PAM["ID"].tolist()
    pam_end = df_PAM["End_pos"].tolist()
    pam_strand = df_PAM["Strand"].tolist()
    pam_coordinate = df_PAM["PAM_pos"].tolist()
    tss_coordinate = [int(i) for i in df_TSS["TSS coordinate"].tolist()[:-1]]

    close_pam_index = []
    for m in pam_index:
        pam_to_search = all_pam[m]
        v = binary_search(tss_coordinate, 0, len(tss_coordinate) - 1, pam_to_search, up_range, down_range)
        if v != -1:
            close_pam_index.append(m)

    short_list = []
    for i in close_pam_index:
        chrom_id = pam_chrom[i]
        pam_start = all_pam[i]
        pam_end_pos = pam_end[i]
        strand = pam_strand[i]
        coordinate = pam_coordinate[i]
        short_list.append([chrom_id, strand, pam_start, pam_end_pos, coordinate])
    df = pd.DataFrame(short_list, columns=['ID', "Strand", "Start_pos", "End_pos", "PAM_pos"])
    print("Found " + str(len(short_list)) + " PAM sequences in proximity to at least one tss" + "\n")

    df.to_csv(shortlist, sep='\t', index=False, header=False)
    print("PAM shortlist saved as: " + shortlist + "\n")
    return shortlist  # returns path of the shortlist of candidates


def sgRNA_finder(sg_candidate, tss, reference, up_range, down_range, sg_outfile):
    global gbk
    df_candidate = pd.read_csv(sg_candidate, sep='\t', names=['ID', "Strand", "Start_pos", "End_pos", "PAM_pos"])
    candidate_start = df_candidate["Start_pos"].tolist()  # start_pos of pam sequences
    candidate_end = df_candidate["End_pos"].tolist()
    candidate_strand = df_candidate["Strand"].tolist()
    candidate_coordinate = df_candidate["PAM_pos"].tolist()
    pam_index = df_candidate.index.tolist()
    t2 = timeit.default_timer()
    df_TSS = pd.read_csv(tss, sep='\t')
    tss_index = df_TSS.index.tolist()[:-1]
    tss_coordinate = [int(i) for i in df_TSS["TSS coordinate"].tolist()[:-1]]
    old_locus = df_TSS["Locus_tag"][:-1].tolist()
    dic_locus = locus_dic(gbk)
    dic_protein = pro_dic(gbk)


    # screen seed regions to avoid off-target effect (12 bp seed regions of the sgRNAs are screened)
    # generates a dictionary {index: sequence}
    with open(reference, "r") as f:
        for record in SeqIO.parse(f, "fasta"):  # read sequence(s) from fasta
            gDNA = str(record.seq)
            candidate_seq = {}
            sg_list = []
            for n in pam_index:
                sequence = str(record.seq[candidate_start[n]:candidate_end[n]])  # sgRNA sequence (default 20 bp)
                seed = sequence[0:12]  # seed region of the sgRNA (12bp)
                rc_sequence = reverse_complement(sequence)
                rc_seed = rc_sequence[0:12]
                if candidate_strand[n] == '+':  # if NGG, keep original sequence
                    if gDNA.count(seed) == 1:
                        if gDNA.count(reverse_complement(seed)) == 1:
                            candidate_seq.update({n: sequence})
                else:  # if CCN, keep reverse complementary
                    if gDNA.count(rc_seed) == 1:
                        if gDNA.count(reverse_complement(rc_seed)) == 1:
                            candidate_seq.update({n: rc_sequence})
            t3 = timeit.default_timer()
            print("Finished screening seed regions in: ", format((t3 - t2), '.2f'), " seconds." + "\n")

            sgRNA = []
            for t in tss_index:
                for m in candidate_seq:
                    if -up_range < candidate_start[m] - tss_coordinate[t] < down_range - 20:
                        strand = candidate_strand[m]
                        old_tag = old_locus[t]
                        new_locus = locus_to_new(old_tag, dic_locus)
                        pro_id = locus_to_protein(old_tag, dic_protein)
                        tss_pos = tss_coordinate[t]
                        sg_start = candidate_start[m]
                        sg_end = candidate_end[m]
                        pam_pos = candidate_coordinate[m]
                        sg_seq = ''.join(['5\'-', candidate_seq[m], '-3\''])
                        sgRNA.append([strand, old_tag, new_locus, pro_id, pam_pos, tss_pos, sg_start, sg_end, sg_seq])

            # Assign ids to each unique sgRNA in sgRNA
            sg_start_list = []
            label = {}
            for entry in sgRNA:
                sg_start_list.append(entry[5])
            temp = defaultdict(lambda: len(temp))
            for ele in sg_start_list:
                label.update({ele: temp[ele]})
            for entry in sgRNA:
                entry.insert(0, ''.join(['SGR_Ab', str(label[entry[5]]+1).zfill(4)]))

            df = pd.DataFrame(sgRNA, columns=['SGR_Ab ID', 'Strand', 'Locus tag (old)', 'Locus tag (new)', 'Protein ID',
                                              'PAM_pos', 'TSS coordinate', 'SGR start', 'SGR end', 'SGR sequence'])
            df.to_csv(sg_outfile, sep='\t', index=False)
    t4 = timeit.default_timer()
    print("Finished recording sgRNA candidates: ", format((t4 - t3), '.2f'), " seconds." + "\n")
    print("Recorded ", len(candidate_seq), "sgRNA sequences." + "\n")
    print("sgRNA sequences saved as: ", sg_outfile, ".", "\n")


def binary_search(arr_tss, l, r, p, up, down):
    if r >= l:
        mid = l + (r - l) // 2;
        if -up <= p - arr_tss[mid] <= down - 20:
            return mid
        elif p - arr_tss[mid] < -up:
            return binary_search(arr_tss, l, mid - 1, p, up, down)
        else:
            return binary_search(arr_tss, mid + 1, r, p, up, down)
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
def pro_dic(annotation_gbk):
    recs = [rec for rec in SeqIO.parse(annotation_gbk, "genbank")]
    dic_protein = {}
    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == "CDS"]
        for feat in feats:
            if 'old_locus_tag' in feat.qualifiers and 'protein_id' in feat.qualifiers:
                locus = ''.join((feat.qualifiers['old_locus_tag']))
                protein_id = ''.join((feat.qualifiers['protein_id']))
            dic_protein.update({locus: protein_id})
        return dic_protein


# Update a given old tag, return a new tag
def locus_to_protein(locus_old, dic_protein):
    if dic_protein.get(locus_old) is None:
        protein_id = ''
    else:
        protein_id = dic_protein.get(locus_old)
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
    gbk = opts.gbk
    up_range = int(opts.up_range)
    down_range = int(opts.down_range)
    sg_outfile = opts.outfile
    base_pam = os.path.splitext(pam_filename)[0]
    shortlist = pam_dir + '/' + 'PAM_shortlist_' + base_pam + '.bed'

    if not (os.path.isabs(sg_outfile)):
        sg_outfile = pam_dir + '/' + sg_outfile

    t0 = timeit.default_timer()
    sg_candidate = tss_pam(pam, tss, up_range, down_range, shortlist)
    t1 = timeit.default_timer()
    print("finished PAM shortlist in: ", format(t1 - t0, '.2f'), " seconds" + "\n")
    sgRNA_finder(sg_candidate, tss, reference, up_range, down_range, sg_outfile)


if __name__ == '__main__':
    main()
