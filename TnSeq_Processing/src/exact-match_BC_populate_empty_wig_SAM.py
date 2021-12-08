## Defne Surujon
## April 5, 2018

# edited by Amy, April 6, 2020
# edited by Yunfei, Dec 07, 2021

# take a wig file for a specific strain, and add the read counts
# from a SAM file for a certain experiment (-S in bowtie)

from Bio import SeqIO
from optparse import OptionParser
from collections import Counter
import re


options = OptionParser(usage='%prog -i [input(empty) wig file] -m [input map] -o [output wig file]',
                       description=
                       """
Specify empty wig file, input map file, and output wig file.
Add the option "--old" to specify the old format for map files
(e.g. output of bowtie 1).
                         """)


options.add_option("-i","--infile",dest="wigfile",
                   help="input .wig file (empty)")
options.add_option("-s","--sam",dest="samfile",
                   help="input .sam file")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output .wig file")

def make_map_dict(samfilename):
    mapdict = Counter()
    f = open(samfilename)
    flines = f.readlines()
    f.close()
    for i in flines:
        if any(header in i for header in ["@HD","@SQ","@PG"]) == False:
            linelist = i.split()
            thisstrand = int(linelist[1]) & 0b00010000
            thispos = int(linelist[3])
            thislen = len(linelist[9])
            if thisstrand == 0: # + strand (flag for reverse complement is off)
                thispos -= 2
            else: # - strand (flag for reverse complement is on)
                thispos += thislen 
            # Yunfei: unlike .map file, .sam file does not provide read count, thus every occurence of read at each site in the empty wig is added up.
            mapdict[thispos] += 1
            #if thispos in mapdict:
            #    mapdict[thispos] = 1
            #else:
            #    mapdict[thispos] += 1
    return(mapdict)

def populate_wig_file(wigfile,infile,outfile):
    # get all TA sites
    f=open(wigfile)
    lines=[line.split() for line in f.readlines()]
    f.close()
    ta_table = {int(line[0]):0 for line in lines[1:]}

    # get occupied TA sited
    mdict = make_map_dict(infile)

    
    print("Total number of TA sites: ", len(ta_table))
    print("Number of occupied TA sites: ",len(mdict))
    # add occupied ta sites to the master list
    # Amy: changed so there are only exact matches
    q = 0
    q3 = 0
    for tasite in mdict:
        try:
            # seems like the positions on the map files were 0-indexed,
            # which is why the +1 correction is there. 
            # Yunfei (2021/12/08): rermoved the +1 correction since SAM files are 1-indexed
            ta_table[tasite] += mdict[tasite]
            q+=1
        except KeyError:
            q3+=1
            
    print("Number of insertions at a TA site: ", q)
    print("Number of insertions not at a TA site: ", q3, "(Discarded)")
    talist = [[key, ta_table[key]] for key in sorted(ta_table.keys())]

    strainname = wigfile[:-4]
    g = open(outfile,'w')
    g.write('variableStep chrom='+strainname+'\n')
    for tasite in talist:
        linestr = [str(tasite[0]), str(tasite[1])]
        g.write("\t".join(linestr)+"\n")
    g.close()



def main():
    opts, args = options.parse_args()
    wigfilename = opts.wigfile
    samfilename = opts.samfile
    outfilename = opts.outputfile

    populate_wig_file(wigfilename,samfilename,outfilename)

if __name__ == '__main__':
    main()
