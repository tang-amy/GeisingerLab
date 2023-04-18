## Yunfei Dai
## 17Apr2023

import os, re, sys

INFILE_DIR = sys.argv[1]
OUTPUT_DIR = sys.argv[2]

for seed in os.listdir(INFILE_DIR):
    if re.match("^seed\d$", seed):
        folder = os.path.join(INFILE_DIR, seed)
        for file in os.listdir(folder):
            if not re.match(r"^.+Input.+$", file):
                exp_name = re.search(r"^BfmR-ChIP-(\d{2})-(\d)(?=.+\.bam$)", file)
                experiment = exp_name.group(1)
                replicate = exp_name.group(2)
                # 17Apr2023 use BfmR-ChIP-28-2 as control
                control = re.search("BfmR-ChIP-28-2"+"Input[A-Za-z0-9_.]+bam", ";".join(list(os.listdir(folder))))
                infile = os.path.join(INFILE_DIR, seed, file)
                control = os.path.join(INFILE_DIR, seed, control.group(0))
                outname = exp_name.group(0) + "." + seed + ".ext_size200"
                outname = os.path.join(OUTPUT_DIR, outname)
                cmd = "macs2 callpeak -t " + \
                     infile + " -c " + control + \
                        " -g 3.8e6 -n " + outname + \
                            " --nomodel --extsize 200"
                #print([os.path.basename(infile), os.path.basename(control)])
                os.system(cmd)
                print("Processing: " + infile)
