## Yunfei Dai
## 05/30/2020


from Bio import SeqIO
import pandas as pd
from bisect import bisect_left


Kroger = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/Kroger_TSS.txt'
Prados = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/Prados_TSS_list.txt'
csv = "/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/ChipSeq-RnaSeq/ab17978-bfmRS_bfmS_vsWT-DGE-resultstolab.csv"
bed_file = "/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/call_peak/macs2_500k_ext200" \
           "/consensus_fasta_bed_intersect/concensus_peak/seed2/BfmR-ChIP-28_seed2.consensus_peak.bed"
gbk = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/NZ_CP012004.gbk'


def read_tss(Kroger, Prados):
    df_TSS_Kroger = pd.read_csv(Kroger, sep='\t', skipfooter=1, engine='python')
    # tss coordinates in both the Korger and Prados tables are unique
    df_TSS_Kroger.loc[df_TSS_Kroger['Primary'] == 1, 'type'] = 'primary'
    df_TSS_Kroger.loc[df_TSS_Kroger['Secondary'] == 1, 'type'] = 'secondary'
    df_TSS_Kroger.loc[df_TSS_Kroger['Internal'] == 1, 'type'] = 'internal'
    df_TSS_Kroger.loc[df_TSS_Kroger['Antisense'] == 1, 'type'] = 'antisense'
    df_TSS_Kroger.loc[df_TSS_Kroger['Orphan'] == 1, 'type'] = 'orphan'
    df_TSS_Prados = pd.read_csv(Prados, sep='\t', engine='python')
    # 'Locus tag' in Prados table is the ORF following each TSS
    return df_TSS_Kroger, df_TSS_Prados


def read_rna_seq(csv):
    df_rna_seq = pd.read_csv(csv)
    df_rna_seq = df_rna_seq[['rank based on p(adj) ∆RSvsWT', 'Gene', 'Geneid', 'Chr', 'Start', 'End',
                             'Strand', 'Length', 'Gene Order', 'Protein ID', 'K Number', 'Global Gene',
                             'Global Product', 'Tag1', 'Tag2', 'Tag3', 'Category1', 'Category2',
                             'Category3', 'Global Path', 'Entry', 'query', 'Entry name', 'Protein names',
                             'Gene names', 'Organism', 'Length.1', 'Gene names  (primary )', 'baseMean',
                             'log2FoldChange-∆RS', 'lfcSE-∆RS', 'stat-∆RS', 'pvalue-∆RS', 'padj-∆RS',
                             'baseMean-∆S', 'log2FoldChange-∆S', 'lfcSE-∆S', 'stat-∆S', 'pvalue-∆S', 'padj-∆S']]
    return df_rna_seq


def get_peak(x, y):
    # find the mid point of a peak
    m = x + (y-x) // 2
    return m


# Correspond old "ACX_60" locus tags to new "ACX_RS60" locus tags using .gbk genome file
# Ignore old tags that cannot be found in the gbk file
# Return a dictionary {old_tag: new_tag}
def locus_dic(annotation_gbk):
    recs = [rec for rec in SeqIO.parse(annotation_gbk, "genbank")]
    dic_locus = {}
    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == "gene"]
        for feat in feats:
            if 'old_locus_tag' in feat.qualifiers:
                old = ''.join((feat.qualifiers['old_locus_tag']))
                new = ''.join((feat.qualifiers['locus_tag']))
            dic_locus.update({old: new})

        dic_locus.update({old: new})
    return dic_locus


def locus_to_new(locus_old, dic_locus):
    # Update a given old tag, return a new tag
    if dic_locus.get(locus_old) is None:
        locus_new = ''
    else:
        locus_new = dic_locus[locus_old]
    return locus_new


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
        if gene_list[match] - coordinate < 500:
            tag = gene_dic_pos[gene_list[match]][2]
        else:
            tag = original_tag
    elif strand == '-':
        gene_list = list(gene_dic_neg.keys())
        match = bisect_left(gene_list, coordinate)
        if coordinate - gene_dic_neg[gene_list[match-1]][1] < 500:
            tag = gene_dic_neg[gene_list[match-1]][2]
        else:
            tag = original_tag
    return tag


def rna_seq_info(locus, df_rna_seq):
    # link information from rna seq table by locus tag
    log2FC_RS = df_rna_seq.loc[df_rna_seq['Gene'] == locus, 'log2FoldChange-∆RS'].item()
    padj_RS = df_rna_seq.loc[df_rna_seq['Gene'] == locus, 'padj-∆RS'].item()
    log2FC_S = df_rna_seq.loc[df_rna_seq['Gene'] == locus, 'log2FoldChange-∆S'].item()
    padj_S = df_rna_seq.loc[df_rna_seq['Gene'] == locus, 'padj-∆S'].item()
    # reformat to scientific notation
    try:
        padj_RS = '{:.2e}'.format(float(padj_RS))
    except ValueError:
        pass
    try:
        padj_S = '{:.2e}'.format(float(padj_S))
    except ValueError:
        pass
    return log2FC_RS, padj_RS, log2FC_S, padj_S


def main():
    df_TSS_Kroger, df_TSS_Prados = read_tss(Kroger, Prados)
    df_rna_seq = read_rna_seq(csv)
    df_peak = pd.read_csv(bed_file, sep='\t', names=['chrom', 'start', 'end', 'intersect', 'files', 's1', 's2', 's3'])
    df_peak.drop(['intersect', 'files', 's1', 's2', 's3'], axis=1, inplace=True)
    df_peak.rename(columns={'start': 'peak_start', 'end': 'peak_end'}, inplace=True)
    dic = locus_dic(gbk)
    gene_dic_pos, gene_dic_neg = gene_dic(gbk)
    result = []
    for peak in df_peak.index.tolist():
        start = df_peak['peak_start'][peak]
        end = df_peak['peak_end'][peak]
        Kroger_tss = df_TSS_Kroger.loc[(df_TSS_Kroger['TSS coordinate'] >= start - 500) &
                                       (df_TSS_Kroger['TSS coordinate'] <= end + 500), 'TSS coordinate'].tolist()

        # Find TSSs from Kroger table
        if not Kroger_tss:
            result.append([peak, start, end, 'Kroger', 'no adjacent Kroger TSS', '', '', '', '', '', '', ''])
        else:
            for match in Kroger_tss:
                strand = df_TSS_Kroger.loc[df_TSS_Kroger['TSS coordinate'] == match, 'Strand'].item()
                tss_type = df_TSS_Kroger.loc[df_TSS_Kroger['TSS coordinate'] == match, 'type'].item()
                original_locus = locus_to_new(df_TSS_Kroger.loc[df_TSS_Kroger['TSS coordinate'] == match,
                                                                'Locus_tag'].item(), dic)
                if tss_type == 'internal':
                    locus = get_tss_gene(match, original_locus, strand, gene_dic_pos, gene_dic_neg)
                elif tss_type == 'antisense':
                    locus = get_tss_gene(match, original_locus, strand, gene_dic_pos, gene_dic_neg)
                    if locus != original_locus:
                        tss_type = 'primary (corrected from antisense)'
                else:
                    locus = original_locus
                if locus in df_rna_seq['Gene'].tolist():
                    log2FC_RS, padj_RS, log2FC_S, padj_S = rna_seq_info(locus, df_rna_seq)
                    result.append([peak, start, end, 'Kroger', match, original_locus, locus, strand, tss_type,
                                   log2FC_RS, padj_RS, log2FC_S, padj_S])
                else:
                    result.append([peak, start, end, 'Kroger', match, original_locus, locus, strand, tss_type,
                                   '', '', '', ''])

        # Find TSSs from Prados table
        Prados_tss = df_TSS_Prados.loc[(df_TSS_Prados['TSS coordinate'] >= start - 500) & (df_TSS_Prados['TSS coordinate'] <= end + 500),
                            'TSS coordinate'].tolist()
        if not Prados_tss:
            result.append([peak, start, end, 'Prados', 'no adjacent Prados TSS', '', '', '', '', '', '', '', ''])
        else:
            for match in Prados_tss:
                original_locus = df_TSS_Prados.loc[df_TSS_Prados['TSS coordinate'] == match, 'Locus tag'].item()
                strand = df_TSS_Prados.loc[df_TSS_Prados['TSS coordinate'] == match, 'Strand'].item()
                if original_locus not in df_rna_seq['Gene'].tolist():
                    locus = get_tss_gene(match, original_locus, strand, gene_dic_pos, gene_dic_neg)
                else:
                    locus = original_locus
                if locus in df_rna_seq['Gene'].tolist():
                    log2FC_RS, padj_RS, log2FC_S, padj_S = rna_seq_info(locus, df_rna_seq)
                    result.append([peak, start, end, 'Prados', match, original_locus, locus, strand, '', log2FC_RS, padj_RS, log2FC_S, padj_S])
                else:
                    result.append([peak, start, end, 'Prados', match, original_locus, locus, strand, '', '', '', '', ''])

    # Write results into output file
    df = pd.DataFrame(result, columns=['peak_index', 'peak_start', 'peak_end', 'source', 'tss_coordinate',
                                       'original_locus_tag', 'Yunfei_adjusted_tag','tss_strand', 'tss_type',
                                       'log2FoldChange-∆RS', 'padj-∆RS', 'log2FoldChange-∆S', 'padj-∆S'])
    df.to_csv('/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/ChipSeq-RnaSeq/test3.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()