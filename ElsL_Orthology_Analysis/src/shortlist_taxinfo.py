## Yunfei Dai
## 2021/10/26

import pandas as pd
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
    with open(infile, 'r') as records:
        title = ['Genome', 'Taxid', 'Organism', 'Superkingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']
        record = []
        out_data = []
        for line in records:
            line = line.strip().split("_")
            genome = line[0]
            taxid = taxid_dict[genome] # search for ncbi taxid from taxid.map
            try:
                handle = Entrez.efetch(db="taxonomy", id=taxid, mode='text', rettype='xml')
                results = Entrez.read(handle)
                for taxon in results:
                    taxid = taxon["TaxId"]
                    rank = taxon["Rank"]
                    name = taxon["ScientificName"]
                    lineage = taxon["Lineage"].split(";")
                    try:
                        superkingdom = lineage[0]
                    except IndexError:
                        superkingdom = ""
                    try:
                        phylum = lineage[1]
                    except IndexError:
                        phylum = ""
                    try:
                        tax_class = lineage[2]
                    except IndexError:
                        tax_class = ""
                    try:
                        order = lineage[3]
                    except IndexError:
                        order = ""
                    try:
                        family = lineage[4]
                    except IndexError:
                        family = ""
                    try:
                        genus = lineage[5]
                    except IndexError:
                        genus = ""
                    record = [genome, taxid, name, superkingdom, phylum, tax_class, order, family, genus]
                    out_data.append(record)
            except Exception:
                pass
        df_output = pd.DataFrame(out_data)
        df_output.to_csv(outfile, header=title, sep='\t', index=False)
            
if __name__ == '__main__':
    main()
