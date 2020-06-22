import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import os
import numpy as np
from sklearn.neighbors import KernelDensity

DIR = '/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/call_peak/macs2_500k_ext200/xls/'
WT = '/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/call_peak/macs2_500k_ext200/xls/' \
        'BfmR-ChIP-28-1_S1_L005_R1_001.s2_.5M.input_28-1.ext200_peaks.xls'
bfmS = '/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/call_peak/macs2_500k_ext200/xls/' \
        'BfmR-ChIP-49-1_S3_L005_R1_001.s2_.5M.input_28-1.ext200_peaks.xls'
D58A = '/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/call_peak/macs2_500k_ext200/xls/' \
        'BfmR-ChIP-29-1_S2_L005_R1_001.s2_.5M.input_28-1.ext200_peaks.xls'

'''
df_28 = pd.DataFrame({'' : []})
for file in os.listdir(DIR):
    path = os.path.join(DIR, file)
    df = pd.read_csv(path, sep='\t', skiprows=range(28), engine='python')
    df['source'] = 'WT'
    if '28' in file:
        df_28.append(df, inplace=True)
print(df_28)
'''

df_WT = pd.read_csv(WT, sep='\t',skiprows=range(28), engine='python')
df_WT['source'] = 'WT'
FE_WT = df_WT['fold_enrichment'].to_numpy()
df_bfmS = pd.read_csv(bfmS, sep='\t',skiprows=range(28), engine='python')
df_bfmS['source'] = 'bfmR'
FE_bfmS = df_bfmS['fold_enrichment'].to_numpy()
df_d58a = pd.read_csv(D58A, sep='\t',skiprows=range(28), engine='python')
df_d58a['source'] = 'D58A'

cdf = pd.concat([df_WT[['source', 'fold_enrichment']], df_bfmS[['source', 'fold_enrichment']],
                 df_d58a[['source', 'fold_enrichment']]])
mdf = pd.melt(cdf, id_vars=['source'], var_name=['FE'])

fig, ax = plt.subplots(2, figsize=[6,8])
sns.boxplot(x=mdf.source, y=mdf.value, ax=ax[0])
ax[0].set_ylabel('Fold Enrichment')
ax[0].set_xlabel('Strain')
ax[0].set_title('Distribution of Fold Enrichment (500k, seed2, input 28-1)')
figname = '/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/call_peak/macs2_500k_ext200/xls/' \
          'FE_distribution_box_s2_.5M_input28-1.eps'
#box.savefig(figname, format='eps')
#sns.boxplot(y='fold_change', hue='FE', data=mdf)
#plt.show()

sns.distplot(df_WT['fold_enrichment'], hist=False, kde=True, kde_kws={'linewidth': 3}, label='WT', ax=ax[1])
sns.distplot(df_bfmS['fold_enrichment'], hist=False, kde=True, kde_kws={'linewidth': 3}, label='bfmS', ax=ax[1])
#sns.distplot(df_d58a['fold_enrichment'], hist=False, kde=True, kde_kws={'linewidth': 3}, label='D58A')
ax[1].set_xlabel('Fold Enrichment')
ax[1].set_ylabel('Density')
ax[1].set_title('Distribution of Fold Enrichment (500k, seed2, input 28-1)')
figname_dist = '/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/call_peak/macs2_500k_ext200/xls/' \
          'FE_distribution_density_s2_.5M_input28-1.eps'
figname_combined = '/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/call_peak/macs2_500k_ext200/xls/' \
          'FE_distribution_s2_.5M_input28-1.eps'
plt.show()
fig.savefig(figname_combined, format='eps')
