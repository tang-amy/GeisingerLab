import sys
import pandas as pd
import numpy as np


def call_effect(condition_code, conditions):
    # based on condition code, make general / specific call
    record = conditions.loc[
        (conditions['RS vs WT'] == condition_code[0]) & 
        (conditions['R vs WT'] == condition_code[1]) & 
        (conditions['S vs WT'] == condition_code[2])
        ]
    general_call = record['general call'].item()
    specific_call = record['specific call'].item()
    effect = [general_call, specific_call]
    return effect

def get_condition_code(significant, fold_change):
    # generate condition code based on fold change and significance
    # yes, pos = 1
    # yes, neg = -1
    # no = 0
    if significant == 'yes':
        if fold_change > 0:
            condition_code = 1
        else: # fold change < 0
            condition_code = -1
    else:
        condition_code = 0
    return condition_code

def main():
    infile = sys.argv[1]
    condition_table = sys.argv[2]
    outfile = sys.argv[3]

    conditions = pd.read_csv(condition_table, sep='\t', index_col=None)
    conditions = conditions.drop_duplicates()

    df_mastertable = pd.read_csv(infile, sep='\t')
    
    # generate new column: condition code
    df_mastertable['condition_code'] = df_mastertable.apply(
        lambda x: (get_condition_code(x['BfmRS_WT_significant'], x['BfmRS_WT_log2(fold_change)']), 
                get_condition_code(x['BfmR_WT_significant'], x['BfmR_WT_log2(fold_change)']), 
                get_condition_code(x['BfmS_WT_significant'], x['BfmS_WT_log2(fold_change)'])
                ),
                axis = 1)
    # generate new column: general call
    df_mastertable['general call'] = df_mastertable.apply(
        lambda x: call_effect(x['condition_code'], conditions)[0],
        axis = 1
    )
    #generate new column: specific call
    df_mastertable['specific call'] = df_mastertable.apply(
        lambda x: call_effect(x['condition_code'], conditions)[1],
        axis = 1
    )
    df_mastertable.to_csv(outfile, sep='\t', index=False)

if __name__ == "__main__":
    main()
