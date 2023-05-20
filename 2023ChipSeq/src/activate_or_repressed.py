import sys
import pandas as pd
import numpy as np


def call_effect(status_code, conditions):
    record = conditions.loc[
        (conditions['RS vs WT'] == status_code[0]) & 
        (conditions['R vs WT'] == status_code[1]) & 
        (conditions['S vs WT'] == status_code[2])
        ]
    general_call = record['general call'].item()
    specific_call = record['specific call'].item()
    effect = [general_call, specific_call]
    return effect

def get_status_code(significant, fold_change):
    if significant == 'yes':
        if fold_change > 0:
            status_code = 1
        else: # fold change < 0
            status_code = -1
    else:
        status_code = 0
    return status_code

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]

    condition_table = "/work/geisingerlab/Yunfei/2023_Apr_ChipSeq_Analysis/RNAseq_Chipseq_mastertable/repressed_or_activated.txt"
    # yes, pos = 1
    # yes, neg = -1
    # no = 0
    conditions = pd.read_csv(condition_table, sep='\t', index_col=None)
    conditions = conditions.drop_duplicates()

    df_mastertable = pd.read_csv(infile, sep='\t')
    df_mastertable['status_code'] = df_mastertable.apply(
        lambda x: (get_status_code(x['BfmRS_WT_significant'], x['BfmRS_WT_log2(fold_change)']), 
                get_status_code(x['BfmR_WT_significant'], x['BfmR_WT_log2(fold_change)']), 
                get_status_code(x['BfmS_WT_significant'], x['BfmS_WT_log2(fold_change)'])
                ),
                axis = 1)

    df_mastertable['general call'] = df_mastertable.apply(
        lambda x: call_effect(x['status_code'], conditions)[0],
        axis = 1
    )

    df_mastertable['specific call'] = df_mastertable.apply(
        lambda x: call_effect(x['status_code'], conditions)[1],
        axis = 1
    )
    df_mastertable.to_csv(outfile, index=False)

if __name__ == "__main__":
    main()
