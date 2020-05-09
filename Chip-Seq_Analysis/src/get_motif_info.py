## Yunfei Dai
## 05/08/2020


"""
From combined.meme output of each meme-chip result, extract predicted motif sequences and their nsites and e value.

Usage:

python3 get_motif.py [DIR] [OUT]

DIR: directory containing all meme-chip results
OUT: directory of output summary file (.tsv format recommended)

"""

import os
import pandas as pd
from sys import argv

DIR = argv[1] if len(argv)>1 else "/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/meme-chip"
OUT = argv[2] if len(argv)>1 else "/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/meme-chip/motif_summary_sorted.tsv"


def get_info(file):
    path = os.path.abspath(file)
    size = path.split('/')[-4][:-12]
    result_name = path.split('/')[-2][:-15] + '_consensus'
    seed = path.split('/')[-2][-16:-15]
    with open(file, 'r') as f:
        line_list = list(enumerate(f, 1))
        INFO = []
        for num, line in line_list:
            if "MOTIF" in line:
                info1 = line_list[num-1][1].split()
                motif = info1[2].split('-')[0]
                method = '-'.join(info1[2].split('-')[1:])
                info2 = line_list[num+1][1].split()
                nsites = info2[7]
                e_val = float(info2[9])
                if e_val != 0.0e+000:
                    INFO.append([size, result_name, seed, motif, method, nsites, e_val])
    df = pd.DataFrame(INFO, columns=['Downsample size', 'File', 'Seed number', 'MOTIF', 'Method', 'nsites', 'E value'])
    df.sort_values(['E value'], ascending=True, inplace=True)
    df['E value'] = df['E value'].map(lambda x: '{:.2e}'.format(x))
    return df


def main():
    for i in os.listdir(DIR):
        full_path = os.path.join(DIR, i)
        if os.path.isdir(full_path):
            path = os.path.join(full_path, 'bed_intersect_consensus_meme_chip')
            for result in os.listdir(path):
                if 'BfmR-ChIP' in result:
                    meme_file = os.path.join(path, result, 'combined.meme')
                    if os.path.isfile(meme_file):
                        # check if combined.meme exists in the folder
                        df = get_info(meme_file)
                        if not os.path.isfile(OUT):
                            df.to_csv(OUT, sep='\t', header='column_names', index=False)
                        else:
                            df.to_csv(OUT, mode='a', sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()