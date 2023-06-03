import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 

infile = "/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/RNAseq_ChipSeq_sorted_20May2023/BfmR-ChIP-49_seed1.master_table.sorted.20May2023.tsv"
outfile = "/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/RNAseq_ChipSeq_sorted_20May2023/BfmR-ChIP-49_seed1.intergenic.distplot.percentbygroup_split.png"

df = pd.read_csv(infile, sep='\t')
df_subset = df.loc[(df['general call'] == 'activated' ) | (df['general call'] == 'repressed') ]
df_subset_intergenic = df_subset.loc[df['match_type'] == 'intergenic']
#df_subset = df.loc[(df['general call'] == 'activated' ) | (df['general call'] == 'repressed') ]c

sns.set_theme(style='whitegrid')

dist_plot = sns.displot(data=df_subset_intergenic, x='distance_to_match', stat='percent', col='general call', common_norm=False, kde=True) # multiple="stack" for stacked plot
dist_plot.set(xlim=(-500, 500))
dist_plot.set_axis_labels("Distance to Nearest ORF")
plt.savefig(outfile, bbox_inches='tight')