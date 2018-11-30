## Amy Tang
## Oct 31 2018

# creates a wig file that contains all of the insertions at each TA site
# the merged wig file can be used for viewing results on IGV

#!/home/bin/python

import os
from collections import Counter
from optparse import OptionParser

from Bio import SeqIO

options = OptionParser(usage='%prog input output',
                       description="Specify input directory and output file")

options.add_option("-i", "--infolder", dest="inDir",
                   help="directory where .wig files being merged located (.wig)")
options.add_option("-o", "--outfile", dest="outfile",
                   help="output file (.wig)")
options.add_option("-g", "--genbank", dest="genbank",
                   help="genbank file for strain name (.gbk)")


def main():
    opts, args = options.parse_args()
    merged_reads = Counter()
    wig_dir = opts.inDir
    outfile = opts.outfile

    # reading the files to create table of all reads per coordinate summed together
    for file in os.listdir(wig_dir):
        if file.endswith(".wig"):
            f = open(wig_dir + file)
            lines = f.readlines()[1:]
            for line in lines:
                coordinate = int(line.split('\t')[0])
                read_count = int(line.split('\t')[1])
                merged_reads[coordinate] += read_count

    # get wig header
    genomefile = opts.genbank
    strain = SeqIO.read(genomefile, "genbank")
    header = 'variableStep chrom=' + strain.name + '\n'

    # outputting the table to output file
    talist = [[key, merged_reads[key]] for key in sorted(merged_reads.keys())]
    g = open(outfile, 'w')
    g.write(header)
    for tasite in talist:
        linestr = [str(tasite[0]), str(tasite[1])]
        g.write("\t".join(linestr) + "\n")
    g.close()


if __name__ == '__main__':
    main()
