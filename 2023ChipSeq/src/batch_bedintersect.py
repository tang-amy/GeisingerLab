## Yunfei Dai
## 18Apr2023

# This script uses bedtools multiinter to find commonintervals among replicates
# Common intervals present in at least n=2 (by default) replicates are extracted for further analysis 

import os, re, sys
import pandas as pd, numpy as np

INPUT_DIR = sys.argv[1]
OUTPUT_DIR = sys.argv[2]

# Use bedtools multiinter to identify common intervals among replicates
# https://bedtools.readthedocs.io/en/latest/content/tools/multiinter.html
def get_replicates(INPUT_DIR, seeds, samples):
    for seed in seeds:
        for sample in samples:
            replicates = []
            for file in os.listdir(INPUT_DIR):
                if sample == 28: # 18Apr2023: for 28, only include replicate 2, and populate the list with 3 x replicate 2
                    if re.match("^.+" + str(sample) + "-2.seed" + str(seed) + ".+$", file):
                        replicates.append(os.path.join(INPUT_DIR, file))
                else:
                    if re.match("^.+" + str(sample) + "-\d.seed" + str(seed) + ".+$", file):
                        replicates.append(os.path.join(INPUT_DIR, file))
            if len(replicates) == 1:
                replicates += replicates
            outfile =os.path.join(OUTPUT_DIR, ("BfmR-ChIP-" + str(sample) + "_seed" + str(seed) + ".intersect.bed"))
            cmd = "bedtools multiinter -cluster -header -i " + " ".join(replicates) + " > " + outfile # 23Feb2023: added -header option so the output contains file path
            os.system(cmd)

def main():
    get_replicates(INPUT_DIR, [1, 4, 7], [28, 49])

if __name__ == "__main__":
    main()