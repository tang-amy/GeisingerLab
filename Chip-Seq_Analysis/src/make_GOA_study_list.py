## Yunfei Dai
## 07/01/2020

"""

This script generates the study list for GO term enrichment analysis.

Usage:
python make_GOA_study_list.py -i [input] -o [output] -p [p value cutoff] -c [column number of p value] -l [column number of locus tag] -k [skip frist n rows] -s [separator (t or ',')]

Options:
-p p value cutoff to selete genes, default is 0.05
-c column number where the p value is (0 based)
-l column number where the gene locus tag is (0 based), default is 1 (second column)
-k how many rows at the begining to skip (header lines) (0 based), default is 0 (skip first row)
-s separator used by the file format (t or ',') 

"""

from optparse import OptionParser
import csv

options = OptionParser()
options.add_option("-i", "--input", dest="infile", help="provide input file")
options.add_option("-o", "--output", dest="outfile", help="provide output file")
options.add_option("-p", "--pvalue", type="float",dest="pvalue", default=0.05, help="p value cutoff")
options.add_option("-c", "--p_col", type="int", dest="p_col", default=39, help="column number of p value")
options.add_option("-l", "--locus_col", type="int",dest="locus_col", default=1, help="column number of locus tag")
options.add_option("-k", "--skip", dest="skip_row", default=0)
options.add_option("-s", "--separator", dest="separator", default='t')

(opts, args) = options.parse_args()
infile = opts.infile
outfile = opts.outfile
p_val = opts.pvalue
p_col = opts.p_col
locus_col = opts.locus_col
skip_row = opts.skip_row
separator = opts.separator

if separator == 't':
    sep = '\t'
elif separator == ',':
    sep = ','

f = open(infile, 'r')
table = csv.reader(f, delimiter=sep)
o = open(outfile, 'w')
for index, line in enumerate(table):
    if index > skip_row:
        if len(line) > p_col:
            if line[p_col] != "NA":
                pval = float(line[p_col])
                locus = line[locus_col]
                if pval <=  p_val and "ACX60_" in locus:
                    o.write(locus + '\n')
f.close()
o.close()
