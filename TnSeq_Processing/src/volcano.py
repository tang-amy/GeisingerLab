
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

df = pd.read_csv('~/Desktop/LDT2-WTATA-resampling.txt', skiprows=range(6), sep='\t', header = 0)
df.drop(['Name', 'Desc', 'Sum Ctrl', 'Sum Exp', 'Sites', 'Mean Ctrl', 'Mean Exp', 'Delta Mean'], axis=1, inplace=True)
df.rename(columns={'#Orf': 'Gene', 'log2FC': 'Fold Change', 'p-value': 'p-value', 'Adj. p-value': 'adjusted p-value'}, inplace=True)
df.set_index('Gene', inplace=True)
df['Inverse_adj_p'] = df.apply(lambda row: 1/ row['adjusted p-value'], axis=1)
print(df['adjusted p-value'].head())
print(df['Inverse_adj_p'].head())
df_ap = pd.DataFrame(df, columns=['Fold Change', 'Inverse_adj_p'])

df_ap.plot(kind='scatter', x='Fold Change', y='Inverse_adj_p')
plt.ylabel('1/p')
plt.xlabel('Fold Change')
plt.show()
#mpl.style.use(['ggplot'])
