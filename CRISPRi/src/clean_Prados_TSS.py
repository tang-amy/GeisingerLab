## Yunfei Dai
## 04/19/2020


"""
Adjust format of the Prados list for easier processing.
"""

import pandas as pd
infile = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/Prados_TSS_raw.csv'
outfile = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/Prados_TSS_list.txt'
df_prados = pd.read_csv(infile)
df_prados.drop(['nuc', 'pval1', 'pval2', 'fisher_pval', 'FDR', 'pattern', 'TATAAT_seq',
               'TATAAT_pos', 'TATAAT_err', 'TTGACA_seq', 'TTGACA_pos', 'TTGACA_err', 'TSS_cluster_size'], 
               axis=1, inplace=True)
df_prados = df_prados[df_prados['is_rTSS'] == 'Yes']
df_prados.reset_index(inplace=True)
index = df_prados.index.tolist()
locus_tag = df_prados['following_ORF_Name'].tolist()
cluster_id = df_prados['TSS_cluster_ID'].tolist()
strand = df_prados['strand'].tolist()
pos = df_prados['position'].tolist()
inside_ORF = df_prados['is_inside_an_ORF'].tolist()
distance_ORF = df_prados['following_ORF_distance'].tolist()
strand_dic = {'Plus':'+', 'Minus':'-'}
ORF_dic = {0:'No', 1:'Yes'}
tss = []
for i in index:
    coordinate = pos[i]
    tag = locus_tag[i]
    s = strand_dic[strand[i]]
    cluster = cluster_id[i]
    distance = distance_ORF[i]
    inside = ORF_dic[inside_ORF[i]]
    tss.append([coordinate, tag, s, cluster, inside, distance])


df = pd.DataFrame(tss, columns=['TSS coordinate', 'Locus tag', 'Strand', 'TSS cluster ID',
                                'is_inside_ORF', 'Distance to ORF'])
df.to_csv(outfile, sep='\t', index=False)

