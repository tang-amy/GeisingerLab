import os
import sys
import pandas as pd

rnaseq_bfmRS = "/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/rnaseq/bfmRS_WT.xlsx"
rnaseq_bfmR = "/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/rnaseq/bfmR_WT.xlsx"
rnaseq_bfmS = "/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/rnaseq/bfmS_WT.xlsx"

def read_rnaseq_table(file):
    df = pd.read_excel(file)
    df = df.set_index("gene_id")
    return df

def main():

    infile = "/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/MACS2_output_500k_ext200_input28-2/peak_stat_21APR2023_distance_corrected/BfmR-ChIP-49_seed1.nearest_orf.tsv"
    outfile = "/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/test_mastertable.tsv"
    
    rnaseq_table = {}
    rnaseq_table.update({"BfmRS_WT" : read_rnaseq_table(rnaseq_bfmRS)})
    rnaseq_table.update({"BfmR_WT" : read_rnaseq_table(rnaseq_bfmR)})
    rnaseq_table.update({"BfmS_WT" : read_rnaseq_table(rnaseq_bfmS)})

    test_table = read_rnaseq_table(rnaseq_bfmRS)

    df_chipseq = pd.read_csv(infile, sep='\t', index_col=0)
    df_chipseq = df_chipseq.set_index("locus_tag")

    col_names = ["log2(fold_change)", "significant", "protein_id", "product"]

    for i,j in rnaseq_table.items():
        df_rnaseq = j
        for row, column in df_chipseq.iterrows():
            for k in col_names:
                new_column_name = i + "_" + k
                try:
                    df_chipseq.loc[row, new_column_name] = df_rnaseq.loc[row,k]
                except:
                    pass
    print(df_chipseq)

if __name__ == "__main__":
        main()
#print(df_chipseq)
