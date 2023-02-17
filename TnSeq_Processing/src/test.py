import pyupset as pyu
import pandas as pd
from pickle import load
output = '/Users/yunfei/circos-0.69-9/crispr_paper_plots/UpSet/test.txt'
infile = open(output, 'rb')
infile.seek(0)
data_dict = load(infile)
infile.close()
pyu.plot(infile)
