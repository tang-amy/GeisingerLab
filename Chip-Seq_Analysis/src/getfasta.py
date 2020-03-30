## Yunfei Dai
## 03/20/2020

# This script generates .fasta files as input for MEME suite, using .narrowPeak files generated from macs2 peak calling.

# The package "pybedtools" is required. If not already installed, execute the following command in terminal:
# pip3 install pybedtools


from optparse import OptionParser
from pybedtools import BedTool


def main():
    options = OptionParser()

    options.add_option("-i", "--input", dest="infile",
                       help="provide input .narrowPeak file")
    options.add_option("-g", "--reference", dest="reference",
                       help="provide genome reference file (in .fasta)")
    options.add_option("-l", "--sequence_length", dest="seq_length", type="int", default=250,
                       help="length of resulting sequence region included for meme analysis")
    options.add_option("-o", "--output", dest="outfile",
                       help="name of the output file")

    (opts, args) = options.parse_args()

    infile = opts.infile
    reference = opts.reference
    seq_length = opts.seq_length//2
    outfile = opts.outfile
    new_bed = str(outfile.replace(".fasta", ".bed"))

    w = open(new_bed, "w+")
    w.write("track name=new bed from narrowPeak\n")
    with open(infile, "r") as f:
        for line in f:
            ls = line.strip().split()
            start = int(ls[1])
            offset = int(ls[9])
            new_start = start + offset - seq_length
            new_end = start + offset + seq_length
            ls[1] = str(new_start)
            ls[2] = str(new_end)
            w.write("\t".join(ls) + "\n")
    w.close()

    mybed = BedTool(new_bed)
    fasta = BedTool(reference)
    mybed = mybed.sequence(fi=fasta, fo=outfile)


if __name__ == '__main__':
    main()
