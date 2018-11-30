## Amy Tang
## November 6 2018

## original code from Defne Surujon (Feb 20 2018)

# make a prot_table file from a genbank
# 9 columns without headers: gene annotation, start, end,
# strand, protein length, nothing, nothing, nothing, locus_tag

#!/home/bin/python

from Bio import SeqIO
from optparse import OptionParser


options = OptionParser(usage='%prog input output',
                       description="Specify input directory and output file")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk)")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.prot_table)")


def main():
    opts, args = options.parse_args()
    infile = opts.inputfile
    outfile = opts.outputfile
    
    mystrain=SeqIO.read(infile,"genbank")
    prot_table = []
    stranddict = {1:'+',-1:'-'}

    for thisfeature in mystrain.features:
        if thisfeature.type=="CDS":
            try:
                thisline = ["-" for i in range(0,9)]
                thisline[0] = thisfeature.qualifiers['product'][0]
                thisline[1] = str(thisfeature.location.start)
                thisline[2] = str(thisfeature.location.end)
                thisline[3] = stranddict[thisfeature.location.strand]
                thisline[4] = int(abs(thisfeature.location.end-thisfeature.location.start)/3)
                thisline[8] = thisfeature.qualifiers['locus_tag'][0]
                
                # removing "<" or ">" if it is occurs in start location coordinate
                if thisline[1].find(">") != -1 or thisline[1].find("<") != -1:
                    fixed = thisline[1].replace(">", "").replace("<", "")
                    thisline[1] = fixed
                # removing "<" or ">" if it is occurs in end location coordinate
                if thisline[2].find(">") != -1 or thisline[2].find("<") != -1:
                    fixed = thisline[2].replace(">", "").replace("<", "")
                    thisline[2] = fixed

                prot_table.append(thisline)
            except KeyError:
                pass

    g = open(outfile,'w')
    for line in prot_table:
        linestr = [str(i) for i in line]
        g.write("\t".join(linestr)+"\n")
    g.close()


if __name__ == '__main__':
    main()
