## Amy Tang
## Nov 1 2018

# for all map files (bowtie1 default output files) within a directory, make wig files for each map file

#!/home/bin/python

import map_to_wig as btw
from optparse import OptionParser
import os
from Bio import SeqIO


options = OptionParser(usage='%prog input output',
                       description="Specify input directory and output file")

options.add_option("-g", "--genome", dest="genomefile",
                   help="genome file (.gbk)")
options.add_option("-i", "--infolder", dest="inDir",
                   help="folder containing .map files")
options.add_option("-o", "--outfolder", dest="outDir",
                   help="folder to write the .wig files to")


def main():
    opts, args = options.parse_args()
    map_dir = opts.inDir
    genomefile = opts.genomefile
    mystrain = SeqIO.read(genomefile, "genbank")
    ta_table = btw.get_ta(mystrain)

    # creating the output file if it does not exist
    if not os.path.isdir(opts.outDir):
        os.mkdir(opts.outDir)
    wig_dir = opts.outDir

    for f in os.listdir(map_dir):
        if f.endswith(".map"):
            outfile = wig_dir + os.path.splitext(os.path.basename(f))[0] + ".wig"
            btw.write_wig(map_dir + f, outfile, ta_table, mystrain)


if __name__ == '__main__':
    main()
