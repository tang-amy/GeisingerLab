## Yunfei Dai
## 04/13/2022

# This script concatenates mutation predictions from breseq results (.html)
# Input is the parent directory that contains multiple folders, with output files for each sample.
# Usage: python cat_breseq.py -i [breseq results] -n [gene to ignore (optional)] -o [output file (optional)]

from os import path, listdir, getcwd
from bs4 import BeautifulSoup
from optparse import OptionParser
from pathlib import Path
import lxml, csv

options = OptionParser()
options.add_option("-i", "--infile", dest="infile", help="provide input directory containing .html results")
options.add_option("-n", "--ignore", dest="ignore", default="///", help="suppress mutations mapped to specified gene")
options.add_option("-o", "--output", dest="outfile", default="", help="specify output file name and directory")

def get_summary(infile, ignore, outfile):
    for item in listdir(infile):
        file=path.join(infile,item)
        if path.isdir(file)==True:
            prediction = path.join(file,"output","index.html")	
            with open(prediction) as fp:
                soup = BeautifulSoup(fp, 'lxml')
                tables = soup.find_all("table")
                table_1 = tables[1]
                headings = [th.get_text() for th in table_1.find_all("th")]

                dataset = []
                for row in table_1.find_all("tr")[2:]:
                    text = [td.get_text() for td in row.find_all("td")]
                    evidence = text[0]
                    seq_id = text[1]
                    position = text[2]
                    mutation = text[3]
                    annotation = text[4]
                    gene = text[5]
                    description = text[6]
                    # Discard mutations mapped to plasmid genes
                    if seq_id == "NZ_CP012004": 
                        if "ACX60_RS00525" not in gene:
                            if ignore not in gene:
                                result = [item, evidence, seq_id, position, mutation, annotation, gene, description]
                                with open(outfile, 'a', newline='') as csvfile:
                                    output_writer = csv.writer(csvfile, delimiter='\t')
                                    output_writer.writerow(result)

def main():
    opts, args = options.parse_args()
    infile = opts.infile
    ignore = opts.ignore
    outfile = opts.outfile
    header = ["Sample", "Evidence", "Seq_ID", "Position", "Mutation", "Annotation", "Gene", "Description"]
    try:
        with open(outfile, 'w', newline='') as csvfile:
            output_writer = csv.writer(csvfile, delimiter='\t')
            output_writer.writerow(header)
    except FileNotFoundError:
        outfile = path.join(infile, "MutationPredictions_all.txt")
        with open(outfile, 'w', newline='') as csvfile:
            output_writer = csv.writer(csvfile, delimiter='\t')
            output_writer.writerow(header)
    get_summary(infile, ignore, outfile)
    print("Output saved as: " + outfile)

if __name__ == '__main__':
    main()
