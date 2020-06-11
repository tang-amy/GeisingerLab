## Yunfei Dai
## 06/09/2020

"""

This script generates the data file for circos visualization from wig files (from transit).

Usage:

python wig_to_circos.py -i [infile]

Optional:

-o [outfile] specify output directory, default is the same directory as the input

-l [linear or log] whether to convert the values to log scale, default is linear

Output example:

#chr start end value [options]
hs5 50 75 0.75

See circos requirements: http://circos.ca/documentation/tutorials/configuration/data_files/

"""

from optparse import OptionParser
import pandas as pd
import numpy as np
import os


options = OptionParser()
options.add_option("-i", "--infile", dest="infile", help="provide input wig file")
options.add_option("-o", "--outfile", dest="outfile", default="default", help="provide output directory")
options.add_option("-l", "--scale", dest='scale', default="linear", help="linear or log")

opts, args = options.parse_args()
wig = os.path.abspath(opts.infile)
scale = opts.scale
basename = os.path.splitext(os.path.basename(wig))[0]
if opts.scale == 'log':
    if opts.outfile == 'default':
        outfile = os.path.split(wig)[0] + '/log_' + basename + '.txt'
    else:
        outfile = opts.outfile
else:
    if opts.outfile == 'default':
        outfile = os.path.split(wig)[0] + '/linear_' + basename + '.txt'
    else:
        outfile = opts.outfile
df = pd.read_csv(wig, sep='\t', names=['start', 'value'], skiprows=1)
df['chrom'] = 'NZ_CP012004'
df['end'] = df['start'].to_numpy() + 1
df = df[['chrom', 'start', 'end', 'value']]
df = df.loc[df['value'] > 0]
if scale == 'log':
    df['value'] = np.log10(df['value'].to_numpy())
    df = df.loc[[df['value'] > 0]]
values = df['value'].tolist()
df.to_csv(outfile, sep=' ', header=True, index=False)
