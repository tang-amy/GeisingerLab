## Amy Tang
## Oct 28 2018

# original code from Defne Surujon
# make a wig file from a genbank and default bowtie1 output file

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
                   help="output file (.prot_table)")


def make_read_counter(bowtiefile):
    reads_counter = Counter()  # dictionary for counting reads, coordinate:count
    # reading the file
    f = open(bowtiefile)
    flines = f.readlines()
    f.close()
    for i in flines:
        linelist = i.split()
        strand = linelist[1]
        coordinate = int(linelist[3])
        length = len(linelist[4])
        if strand == "-":
            coordinate += length - 2
        reads_counter[coordinate] += 1
    return dict(reads_counter)

def get_ta(strain):
    seq = strain.seq
    seqstr = str(seq).upper()
    taF = [i.start() + 1 for i in re.finditer("TA", seqstr)]
    ta_ALL = taF
    ta_ALL.sort()
    return {i: 0 for i in ta_ALL}


def main():
    opts, args = options.parse_args()
    genomefile = opts.genomefile
    infile = opts.inputfile
    outfile = opts.outputfile
    # get all TA sites
    mystrain = SeqIO.read(genomefile, "genbank")
    ta_table = get_ta(mystrain)

    # get occupied TA sited
    mdict = make_read_counter(infile)
    print(len(ta_table), len(mdict))
    # add occupied ta sites to the master list
    q = 0
    q1 = 0
    q2 = 0
    for tasite in mdict:
        try:
            # seems like the positions on the map files were 0-indexed,
            # which is why the +1 correction is there.
            ta_table[tasite + 1] += mdict[tasite]
        except KeyError:
            if tasite in ta_table.keys():
                ta_table[tasite] += mdict[tasite]
                q1 += 1
            elif tasite + 2 in ta_table.keys():
                ta_table[tasite + 2] += mdict[tasite]
                q1 += 1
            elif tasite - 1 in ta_table.keys():
                ta_table[tasite - 1] += mdict[tasite]
                q2 += 1
            elif tasite + 3 in ta_table.keys():
                ta_table[tasite + 3] += mdict[tasite]
                q2 += 1
            else:
                # ta_table[tasite+1] = mdict[tasite]
                q += 1
            # print("WARNING: Site "+str(tasite)+"does not appear in the genome!!")
    print(len(ta_table), q, q1, q2)

    talist = [[key, ta_table[key]] for key in sorted(ta_table.keys())]

    g = open(outfile, 'w')
    g.write('variableStep chrom=' + mystrain.name + '\n')
    for tasite in talist:
        linestr = [str(tasite[0]), str(tasite[1])]
        g.write("\t".join(linestr) + "\n")
    g.close()


if __name__ == '__main__':
    main()