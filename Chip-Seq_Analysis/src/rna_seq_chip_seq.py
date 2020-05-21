
from Bio import SeqIO
import pandas as pd
from bisect import bisect_left


tss = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/Kroger_TSS.txt'
csv = "/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/ChipSeq-RnaSeq/ab17978-bfmRS_bfmS_vsWT-DGE-resultstolab.csv"
bed_file = "/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/call_peak/macs2_500k_ext200" \
           "/consensus_fasta_bed_intersect/concensus_peak/seed2/BfmR-ChIP-28_seed2.consensus_peak.bed"
gbk = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/NZ_CP012004.gbk'


def read_tss(tss):
    df_TSS = pd.read_csv(tss, sep='\t', skipfooter=1, engine='python')
    df_TSS.set_index('TSS coordinate', inplace=True)
    # tss coordinates in the Korger table are unique
    df_TSS.loc[df_TSS['Primary'] == 1, 'type'] = 'primary'
    df_TSS.loc[df_TSS['Secondary'] == 1, 'type'] = 'secondary'
    df_TSS.loc[df_TSS['Internal'] == 1, 'type'] = 'internal'
    df_TSS.loc[df_TSS['Antisense'] == 1, 'type'] = 'antisense'
    df_TSS.loc[df_TSS['Orphan'] == 1, 'type'] = 'orphan'
    return df_TSS


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
    return dic_locus


def locus_to_new(locus_old, dic_locus):
    # Update a given old tag, return a new tag
    if dic_locus.get(locus_old) is None:
        locus_new = ''
    else:
        locus_new = dic_locus[locus_old]
    return locus_new


def find_tss_match(i, tss_pos):
    # find the closet TSS for a given peak position
    if i in tss_pos:
        match = i
    else:
        insert_pos = bisect_left(tss_pos, i)
        if i - tss_pos[insert_pos-1] <= tss_pos[insert_pos] - i:
            match = tss_pos[insert_pos - 1]
        else:
            match = tss_pos[insert_pos]
    return match


def main():
    df_TSS = read_tss(tss)
    df_rna_seq = read_rna_seq(csv)
    df_peak = pd.read_csv(bed_file, sep='\t', names=['chrom', 'start', 'end', 'intersect', 'files', 's1', 's2', 's3'])
    df_peak.drop(['intersect', 'files', 's1', 's2', 's3'], axis=1, inplace=True)
    df_peak['nearest_tss'] = df_peak.apply(lambda x: find_tss_match((x.start+x.end)//2, df_TSS.index.tolist()), axis=1)
    df_peak.rename(columns={'start': 'peak_start', 'end': 'peak_end'}, inplace=True)
    df_peak['locus_tag'] = [df_TSS['Locus_tag'][i] for i in df_peak['nearest_tss']]
    df_peak['tss_type'] = [df_TSS['type'][i] for i in df_peak['nearest_tss']]
    dic_loc = locus_dic(gbk)
    for i in df_peak.index:
        tag = locus_to_new(df_peak['locus_tag'][i], dic_loc)
        if tag in df_rna_seq['Gene'].tolist():
            df_peak.loc[i, 'log2FoldChange-∆RS'] = df_rna_seq.loc[df_rna_seq['Gene'] == tag, 'log2FoldChange-∆RS'].item()
            df_peak.loc[i, 'pvalue-∆RS'] = df_rna_seq.loc[df_rna_seq['Gene'] == tag, 'pvalue-∆RS'].item()
            df_peak.loc[i, 'log2FoldChange-∆S'] = df_rna_seq.loc[df_rna_seq['Gene'] == tag, 'log2FoldChange-∆S'].item()
            df_peak.loc[i, 'pvalue-∆S'] = df_rna_seq.loc[df_rna_seq['Gene'] == tag, 'pvalue-∆S'].item()


    print(df_peak.head())
    print(df_peak.columns)
    #print(df_TSS.columns)
    #print(df_rna_seq.columns)


if __name__ == '__main__':
    main()