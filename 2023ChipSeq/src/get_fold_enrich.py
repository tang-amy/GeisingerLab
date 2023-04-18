import sys, os
import pandas as pd

def get_fold_enrich(gene, locus, infile):
    fname = os.path.basename(infile)
    fname = fname.split(".")[0]
    df = pd.read_csv(infile, sep='\t')
    try:
        df_match = df.loc[df["locus_tag"] == locus]
        for index, row in df_match.iterrows():
            start = row["start"]
            end = row["end"]
            match_type = row["match_type"]
            fold_enrich = row["average_fold_enrichment"]
            print([gene, locus, fold_enrich])
    except:
        print("no match for " + locus)

def main():
    genes = sys.argv[1]
    indir = sys.argv[2]
    outfile = ""
    for fname in os.listdir(indir):
        if ".tsv" in fname:
            print(fname)
            infile = os.path.join(indir, fname)
            with open(genes, 'r') as f:
                for line in f:
                    line = line.strip()
                    gene = line.split("\t")[0]
                    locus = line.split("\t")[1]
                    get_fold_enrich(gene, locus, infile)

if __name__ == "__main__":
    main()