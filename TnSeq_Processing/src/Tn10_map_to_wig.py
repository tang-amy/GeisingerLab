# Amy Tang
# April 12, 2020

## based on exact-match_BC_populate_empty_wig.py
## original script from Defne Surujon (April 5, 2018)

# add the read counts from a map file for a certain experiment
# and output the read counts for insertion sites in a wig format

from optparse import OptionParser
import re


options = OptionParser(usage='%prog -i [input(empty) wig file] -m [input map] -o [output wig file]',
                       description=
                       """
Specify input map file, output wig file, and strain name.
Add the option "--old" to specify the old format for map files
(e.g. output of bowtie 1).
                         """)

options.add_option("-m","--map",dest="mapfile",
                   help="input .map file")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output .wig file")
options.add_option("-s", dest="strain", 
                   help="strain name")
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
        # adjustment for - strand
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
        # adjustment for - strand
        if thisstrand == "-":
            thispos = thispos + thislen
        
        if thispos in mapdict:
            mapdict[thispos] += thiscount
        else:
            mapdict[thispos] = thiscount
    return(mapdict)

def write_wig_file(strainname, infile, outfile, isold):
  # get occupied sites
    if isold:
        mdict = make_oldmap_dict(infile)
    else:
        mdict = make_map_dict(infile)
    
    # sort the insertion counts by position into a 2D array (2 X unique insertion sites) 
    # containing insertion site, # of insertions at that site
    insertion_counts = [[key, mdict[key]] for key in sorted(mdict.keys())]
    
    # write into the output file
    f = open(outfile, 'w')
    for site in insertion_counts:
      linestr = [str(site[0]), str(site[1])]
      f.write("\t".join(linestr)+"\n")
    f.close()

def main():
    opts, args = options.parse_args()
    mapfilename = opts.mapfile
    outfilename = opts.outputfile
    isold = opts.oldmap
    strainname = opts.strain
    write_wig_file(strainname, mapfilename, outfilename, isold)



if __name__ == '__main__':
    main()
