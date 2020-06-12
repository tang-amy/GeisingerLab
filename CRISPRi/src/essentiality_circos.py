## Yunfei Dai
## 06/12/2020

"""
This script generates data files for circos visualization of essential genes.

1. Finds gene coordinates based on locus tag in the table, using gbk file;

2. If essential, value = 1; Otherwise value = 0;

3. Write chromosome name / start / end / value into output .txt file.
"""

import pandas as pd
from Bio import SeqIO

file = '/Users/yunfei/circos-0.69-9/crispr_paper_circos/TableS1-TnSeq-essentiality.csv'
gbk = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/NZ_CP012004.gbk'
output_mariner = '/Users/yunfei/circos-0.69-9/crispr_paper_circos/' \
                 'data_files_for_plotting/circos_mariner_essentiality.txt'
output_tn10 = '/Users/yunfei/circos-0.69-9/crispr_paper_circos/' \
              'data_files_for_plotting/circos_tn10_essentiality.txt'
col_names = ['locus', 'gene', 'protein ID', 'genbank prot description', 'k',
             'n', 'r', 's', 'zbar', 'mariner_call', 'k.1', 'n.1', 'r.1', 'ovr', 'lenovr',
             'pval', 'padj', 'tn10_call', 'both_call', 'old_A1S_locus',
             'old A1S_Essentiality (Wang et al 2014)', 'AB5075-UW ortholog',
             'AB5075-UW Essentiality (Gallagher et al 2015)']

recs = [rec for rec in SeqIO.parse(gbk, 'genbank')]
gbk_list = {}
for rec in recs:
    feats = [feat for feat in rec.features if feat.type == 'gene']
    for feat in feats:
        tag = ''.join((feat.qualifiers['locus_tag']))
        start = int(feat.location.start)
        end = int(feat.location.end)
        gbk_list.update({tag: (start, end)})

df_input = pd.read_csv(file, sep=',', skiprows=range(2), names=col_names, skipfooter=27, engine='python')
df_input['chrom'] = 'NZ_CP012004'
df_input['start'] = [gbk_list[i][0] for i in df_input['locus'].tolist()]
df_input['end'] = [gbk_list[i][1] for i in df_input['locus'].tolist()]
df_input.loc[df_input['mariner_call'] != 'E', 'mariner_call'] = 0
df_input.loc[df_input['mariner_call'] == 'E', 'mariner_call'] = 1
df_input.loc[df_input['tn10_call'] != 'Essential', 'tn10_call'] = 0
df_input.loc[df_input['tn10_call'] == 'Essential', 'tn10_call'] = 1
df_mariner = df_input[['chrom', 'start', 'end', 'mariner_call']]
df_tn10 = df_input[['chrom', 'start', 'end', 'tn10_call']]

df_mariner.to_csv(output_mariner, sep=' ', header=True, index=False)
df_tn10.to_csv(output_tn10, sep=' ', header=True, index=False)
