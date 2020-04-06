## Defne Surujon
## April 5, 2018

# edited by Amy, April 6, 2020

# take a wig file for a specific strain, and add the read counts
# from a map file for a certain experiment

from Bio import SeqIO
from optparse import OptionParser
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
options.add_option("-m","--map",dest="mapfile",
                   help="input .map file")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output .wig file")
options.add_option("--old", dest = "oldmap",
                   action="store_true", default = False)

def make_map_dict(mapfilename):
    mapdict = {}
    f = open(mapfilename)
    flines = f.readlines()
    f.close()
    for i in flines:
        linelist = i.split()
        thiscount = int(linelist[0])
        thisstrand = linelist[1]
        thispos = int(linelist[2])
        thislen = int(linelist[3])
        // Amy: changed how position is determined
        if thisstrand == "+":
            thispos = thispos - 2
        if thisstrand == "-":
            thispos = thispos + thislen
        
        if thispos in mapdict:
            mapdict[thispos] += thiscount
        else:
            mapdict[thispos] = thiscount
    return(mapdict)

def make_oldmap_dict(mapfilename):
    mapdict = {}
    f = open(mapfilename)
    flines = f.readlines()
    f.close()
    for i in flines:
        linelist = i.split()
        thiscount = int(linelist[0].split('-')[1])
        thisstrand = linelist[1]
        thispos = int(linelist[3])
        thislen = len(linelist[5])
        // Amy: changed how position is determined
        if thisstrand == "+":
            thispos = thispos - 2
        if thisstrand == "-":
            thispos = thispos + thislen
        
        if thispos in mapdict:
            mapdict[thispos] += thiscount
        else:
            mapdict[thispos] = thiscount
    return(mapdict)

def populate_wig_file(wigfile,infile,outfile,isold):
    # get all TA sites
    f=open(wigfile)
    lines=[line.split() for line in f.readlines()]
    f.close()
    ta_table = {int(line[0]):0 for line in lines[1:]}

    # get occupied TA sited
    if isold:
        mdict = make_oldmap_dict(infile)
    else:
        mdict = make_map_dict(infile)

    
    print("Total number of TA sites: ", len(ta_table))
    print("Number of occupied TA sites: ",len(mdict))
    # add occupied ta sites to the master list
    q = 0
    q1 = 0
    q2 = 0
    q3 = 0
    for tasite in mdict:
        try:
            # seems like the positions on the map files were 0-indexed,
            # which is why the +1 correction is there. 
            ta_table[tasite+1] += mdict[tasite]
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
    mapfilename = opts.mapfile
    outfilename = opts.outputfile
    isold = opts.oldmap

    populate_wig_file(wigfilename,mapfilename,outfilename,isold)


if __name__ == '__main__':
    main()
