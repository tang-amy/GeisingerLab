## Yunfei Dai
## 2021/10/26

from Bio import Entrez
from optparse import OptionParser
Entrez.email = "dai.yun@northeastern.edu"

options = OptionParser()
options.add_option("-i", "--infile", dest="infile", help="input is a list of GT genome identifiers")
options.add_option("-r", "--reference", dest="reference", default="/work/geisingerlab/Yunfei/WebOfLife/taxonomy/taxid.map", help="reference is taxid.map")
options.add_option("-o", "--outfile", dest="outfile", help="specify the output file directory")


def get_tax_dict(tax_map):
    # create a dictionary of {genome : taxid} from taxid.map
    with open(tax_map, 'r') as genome_map:
        taxid_dict = {}
        for line in genome_map:
            line = line.strip()
            line = line.split()
            genome = line[0]
            taxid = line[1]
            taxid_dict.update({genome: taxid})
    return taxid_dict


def main():
    opts, args = options.parse_args()
    infile = opts.infile
    outfile = opts.outfile
    reference = opts.reference
    
    taxid_dict = get_tax_dict(reference)
    writer = open(outfile, 'a+')
    with open(infile, 'r') as records:
        for line in records:
            line = line.strip().split("_")
            genome = line[0]
            taxid = taxid_dict[genome] # search for ncbi taxid from taxid.map
            handle = Entrez.efetch(db="taxonomy", id=taxid, mode='text', rettype='xml')
            results = Entrez.read(handle)
            for taxon in results:
                taxid = taxon["TaxId"]
                name = taxon["ScientificName"]
                line_write = '\t'.join([genome, taxid, name])
                writer.write(line_write + '\n')
            
if __name__ == '__main__':
    main()
