import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 

infile = sys.argv[1]
outfile = sys.argv[2]

df = pd.read_csv(infile, sep='\t')
df_subset = df.loc[(df['general call'] == 'activated' ) | (df['general call'] == 'repressed') ]
df_subset_intergenic = df_subset.loc[df['match_type'] == 'intergenic']
#df_subset = df.loc[(df['general call'] == 'activated' ) | (df['general call'] == 'repressed') ]c

sns.set_theme(style='whitegrid')

dist_plot = sns.displot(data=df_subset_intergenic, x='distance_to_match', stat='percent', col='general call', common_norm=False, kde=True) # multiple="stack" for stacked plot
dist_plot.set(xlim=(-500, 500))
dist_plot.set_axis_labels("Distance to Nearest ORF")
plt.savefig(outfile, bbox_inches='tight')
