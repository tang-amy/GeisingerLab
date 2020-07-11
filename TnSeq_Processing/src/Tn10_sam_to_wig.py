# Amy Tang
# July 11, 2020

## based on exact-match_BC_populate_empty_wig.py
## original script from Defne Surujon (April 5, 2018)

# add the read counts from a SAM file for a certain experiment
# and output the read counts for insertion sites in a wig format

from optparse import OptionParser
import re


options = OptionParser(usage='%prog -i [input(empty) wig file] -m [input map] -o [output wig file]',
                       description=
                       """
Specify input sam file, output wig file, and strain name.
                         """)

options.add_option("-S","--sam",dest="samfile",
                   help="input .map file")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output .wig file")
options.add_option("-s", dest="strain", default='NZ_CP012004', 
                   help="strain name")

def make_sam_dict(samfilename):
    samdict = {}
    f = open(samfilename)
    flines = f.readlines()
    f.close()
    for i in flines:
        linelist = i.split()
        strand = int(linelist[2]) & 0b00010000
        if strand == 0:
            thisstrand = "+" # + strand << flag for reverse complement is off
        else:
            thisstrand = "-" # - strand << flag for reverse complement is on
        thiscoordinate = int(linelist[4])
        thislen = len(linelist[10])
        # adjustment for - strand
        if thisstrand == "-":
            thispos = thispos + thislen
        # adding this alignment to the dictionary
        if thispos in samdict:
            samdict[thispos] += 1
        else:
            samdict[thispos] = 1
    return(samdict)

def write_wig_file(strainname, infile, outfile):
    # get occupied sites
    mdict = make_sam_dict(infile)
    
    # create dict of all the positions and set the counts to 0
    insertions_dict = {}
    for i in range(1,3857743+1):
      insertions_dict[i] = 0
    
    print("Number of unique insertion sites:", len(mdict))
    # adding in the counts from the sam file to the insertion dict
    for site in mdict:
      insertions_dict[site] += mdict[site]
    
    # write into the output file
    f = open(outfile, 'w')
    f.write('variableStep chrom='+strainname+'\n')
    # go through all insertion sites and append site, read count into wig output file
    for key in sorted(insertions_dict.keys()):
      linestr = str(key) + "\t" + str(insertions_dict[key])
      f.write(linestr+"\n")
    f.close()

def main():
    opts, args = options.parse_args()
    samfilename = opts.samfile
    outfilename = opts.outputfile
    strainname = opts.strain
    write_wig_file(strainname, samfilename, outfilename)



if __name__ == '__main__':
    main()
