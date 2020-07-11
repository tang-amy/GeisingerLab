## Amy Tang
## Apr 4 2020

# original code from Defne Surujon
# make a wig file from a genbank and SAM file

#!/home/bin/python

from Bio import SeqIO
from optparse import OptionParser
import re
from collections import Counter

options = OptionParser(usage='%prog input output',
                       description="Specify input directory and output file")

options.add_option("-g", "--genome", dest="genomefile",
                   help="genome file (.gbk)")
options.add_option("-i", "--infile", dest="inputfile",
                   help="input map file (.map)")
options.add_option("-o", "--outfile", dest="outputfile",
                   help="output file (.wig)")


def make_read_counter(samfile):
    reads_counter = Counter()  # dictionary for counting reads, coordinate:count
    # reading the file
    f = open(samfile)
    flines = f.readlines()
    f.close()
    for i in flines:
        linelist = i.split()
        strand = int(linelist[2]) & 0b00010000
        coordinate = int(linelist[4])
        length = len(linelist[10])
        if strand == 0: # + strand << flag for reverse complement is off
            coordinate += length
        else: # - strand << flag for reverse complement is on
            coordinate -= 2 # shift the coordinate to the beginning of the TA site
        reads_counter[coordinate] += 1
    return dict(reads_counter)

def get_ta(strain):
    seq = strain.seq
    seqstr = str(seq).upper()
    taF = [i.start() + 1 for i in re.finditer("TA", seqstr)]
    ta_ALL = taF
    ta_ALL.sort()
    return {i: 0 for i in ta_ALL}


def write_wig(f, outfile, ta_table, mystrain):
    # get occupied TA sited
    mdict = make_read_counter(f)
    print(len(ta_table), " total TA sites found\n")
    print(len(mdict), "inserted TA sites found\n")

    # add occupied ta sites to the master list
    q = 0 # not matched within error range (+-2)
    q1 = 0 # matched within +- 1 error
    q2 = 0 # matched within +- 2 error
    for tasite in mdict:
        try:
            # SAM files are 1-indexed so no offset is necessary
            ta_table[tasite] += mdict[tasite]
        except KeyError:
            if tasite - 1 in ta_table.keys():
                ta_table[tasite - 1] += mdict[tasite]
                q1 += 1
            elif tasite + 1 in ta_table.keys():
                ta_table[tasite + 1] += mdict[tasite]
                q1 += 1
            elif tasite - 2 in ta_table.keys():
                ta_table[tasite - 2] += mdict[tasite]
                q2 += 1
            elif tasite + 2 in ta_table.keys():
                ta_table[tasite + 2] += mdict[tasite]
                q2 += 1
            else:
                q += 1
    print(q, " no match within +/- 2 bp\n")
    print(q1, " mapped within +/- 1 bp\n")
    print(q2, " mapped within +/- 2 bp\n")
    print("Finished")

    talist = [[key, ta_table[key]] for key in sorted(ta_table.keys())]

    g = open(outfile, 'w')
    g.write('variableStep chrom=' + mystrain.name + '\n')
    for tasite in talist:
        linestr = [str(tasite[0]), str(tasite[1])]
        g.write("\t".join(linestr) + "\n")
    g.close()


def main():
    opts, args = options.parse_args()
    genomefile = opts.genomefile
    infile = opts.inputfile
    outfile = opts.outputfile
    # get all TA sites
    mystrain = SeqIO.read(genomefile, "genbank")
    ta_table = get_ta(mystrain)

    write_wig(infile, outfile, ta_table, mystrain)


if __name__ == '__main__':
    main()
