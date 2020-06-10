
from sys import argv
from Bio import SeqIO
from bisect import bisect_left

gbk = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/NZ_CP012004.gbk'

def gene_dic(gbk):
    # Parse gbk file to get dictionaries containing gene positions and locus tags
    recs = [rec for rec in SeqIO.parse(gbk, "genbank")]
    gene_dic_pos = {}
    gene_dic_neg = {}
    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == "gene"]
        for feat in feats:
            if feat.location.strand == 1:
                pos = int(feat.location.start)
                end = int(feat.location.end)
                tag = ''.join((feat.qualifiers['locus_tag']))
                gene_dic_pos.update({pos: [pos, end, tag]})
            elif feat.location.strand == -1:
                pos = int(feat.location.start)
                end = int(feat.location.end)
                tag = ''.join((feat.qualifiers['locus_tag']))
                gene_dic_neg.update({pos: [pos, end, tag]})
    return gene_dic_pos, gene_dic_neg


def get_tss_gene(coordinate, original_tag, strand, gene_dic_pos, gene_dic_neg):
    # Return the locus tag of the gene whose transcription is regulated by a given TSS
    if strand == '+':
        gene_list = list(gene_dic_pos.keys())
        match = bisect_left(gene_list, coordinate)
        print(gene_dic_pos[gene_list[match]])
        if gene_list[match] - coordinate < 500:
            tag = gene_dic_pos[gene_list[match]][2]
        else:
            tag = original_tag
    elif strand == '-':
        gene_list = list(gene_dic_neg.keys())
        match = bisect_left(gene_list, coordinate)
        if coordinate - gene_dic_neg[gene_list[match-1]][1] < 500:
            tag = gene_dic_neg[gene_list[match-1]][2] + '_new'
        else:
            tag = original_tag
    return tag


gene_dic_pos, gene_dic_neg = gene_dic(gbk)
tag1 = get_tss_gene(3646, 'ori', '+', gene_dic_pos, gene_dic_neg)
tag2 = get_tss_gene(104781, 'ori', '-', gene_dic_pos, gene_dic_neg)
print('tag1: ', tag1)
print('tag2: ', tag2)