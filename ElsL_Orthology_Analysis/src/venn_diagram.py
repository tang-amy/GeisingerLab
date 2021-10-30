## Yunfei Dai
## 10/08/2021

from sys import argv
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import pandas as pd

infile_1 = argv[1]
infile_2 = argv[2]
infile_3 = argv[3]
outfile = argv[4]
try:
    keyword = argv[5]
except Exception:
    keyword = "Organism"

df_1 = pd.read_csv(infile_1, sep='\t')
df_2 = pd.read_csv(infile_2, sep='\t')
df_3 = pd.read_csv(infile_3, sep='\t')

subset_1 = df_1[keyword].tolist()
subset_2 = df_2[keyword].tolist()
subset_3 = df_3[keyword].tolist()

Abc = len(set(subset_1))
aBc = len(set(subset_2))
abC = len(set(subset_3))

ABc = len(list(set(subset_1).intersection(subset_2)))
AbC = len(list(set(subset_1).intersection(subset_3)))
aBC = len(list(set(subset_2).intersection(subset_3)))

ABC = len(list(set(subset_1).intersection(subset_2, subset_3)))

venn3(subsets = (Abc, aBc, ABc, abC, AbC, aBC, ABC), set_labels = ("ElsL", "LdcA", "LdcV"), alpha=0.5)
plt.savefig(outfile, format='pdf')

print(infile_1 + "and" + infile_2 + ": "  + str(list(set(subset_1).intersection(subset_2))))
print(infile_1 + "and" + infile_3 + ": "  + str(list(set(subset_1).intersection(subset_3))))
print(infile_2 + "and" + infile_3 + ": "  + str(list(set(subset_2).intersection(subset_3))))
