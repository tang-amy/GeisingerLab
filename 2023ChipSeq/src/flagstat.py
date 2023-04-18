## Yunfei Dai
## 27NOV2022

# This script counts the number of alignments mapped to each reference and saves the results to a tsv file

import os, pandas as pd
import subprocess
from collections import Counter
SAM_DIR = "/work/geisingerlab/Yunfei/2022_Nov_ChipSeq_Analysis/Mapped_SAM"
outfile = "/work/geisingerlab/Yunfei/2022_Nov_ChipSeq_Analysis/flagstat.tsv"

def flag_stat(infile):
    fname = os.path.basename(infile).split(".")[0]
    samfile = subprocess.Popen(["samtools", "view", infile], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, bufsize=1)
    count1 = 0 # NZ_CP012004.1 (chromosome)
    count2 = 0 # pAB1:  NC_009083.1
    count3 = 0 # pAB2:  NC_009084.1
    count4 = 0 # pAB3:  NZ_CP012005.1
    mapping = []
    for line in samfile.stdout:
        dline = line.decode('ascii')
        if not dline.startswith('@'):
            reference = dline.split()[2]
            mapping.append(reference)
    occurrence = Counter(mapping)
    count1 = occurrence["NZ_CP012004.1"] # chromosome
    count2 = occurrence["NC_009083.1"] # pAB1
    count3 = occurrence["NC_009084.1"] # pAB2
    count4 = occurrence["NZ_CP012005.1"] # pAB3
    unmapped = occurrence["*"] # unmapped
    total = sum([count1, count2, count3, count4, unmapped])
    total_mapped = sum([count1, count2, count3, count4])
    return [fname, total, total_mapped, count1, count2, count3, count4, unmapped]

header = ["File", "Total reads", "Total mapped", "NZ_CP012004.1", "NC_009083.1", "NC_009084.1", "NZ_CP012005.1", "Unmapped"]
writer = []

for file in os.listdir(SAM_DIR):
    infile = os.path.join(SAM_DIR, file)
    writer.append(flag_stat(infile))
df = pd.DataFrame(writer, columns=header)
print(df.head())
#df.to_csv(outfile, sep='\t', index=False)