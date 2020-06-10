## 06/09/2020

"""
1) sum all of the numerical values in columns 2 through 24 (This will give me the sum of all reads across all of the Tn10 subpools)

2) divide by a few different numbers: divide by 500, OR by 1000, OR by 1500, by 2000, OR by 2500. (This will get the normalized read values into a range that is better for TRANSIT -- please output the result as a separate file for each, as you will see in step 3)

3) for each, output as a 2-column text file with just column 1 and the value resulting after the division. add the header
"variableStep chrom=NZ_CP012004" so it can be used in transit and igv as wig

Numpy and pandas (version 0.24.0 or above) are required.

Usage:

python normalize_tn10_combined_wig_files.py [INPUT] [OUTPUT DIRECTORY]

"""
import pandas as pd
import numpy as np
from sys import argv
import os

INPUT = argv[1] if len(argv)>1 else "/Volumes/Seagate/TN10-TTR-combined_wig_files/TN10-TTR-combined_wig_output-20200609.txt"
OUTPUT = argv[2] if len(argv)>2 else "/Volumes/Seagate/TN10-TTR-combined_wig_files"
base = '.'.join(os.path.abspath(INPUT).split('/')[-1].split('.')[:-1])

df = pd.read_csv(INPUT, sep='\t', skiprows=[i for i in range(0, 27)],
                 engine='python', header=None)
v = df.iloc[:, 1:24].to_numpy(dtype='float')
total = np.sum(v, axis=1)
for i in [500, 1000, 1500, 2000, 2500]:
    outname = 'denominator_' + str(i) + '_' + base + '.wig'
    output = os.path.join(OUTPUT, outname)
    with open(output, 'w+') as f:
        f.write('variableStep chrom=NZ_CP012004\n')
    division = np.true_divide(total, i)
    df_result = pd.DataFrame({'col1': df.iloc[:, 0], 'col2': division})
    df_result.to_csv(output, sep='\t', index=False, header=False, mode='a')

