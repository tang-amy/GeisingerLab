#!/bin/bash

## Yunfei Dai
## 04/28/2020

<<Comment
  Run macs2 with downsample treatment files against input.
  Usage:
  bash batch_call_peak.sh <input directory> <output directory>
Comment
SRC="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DIR=${1?Error: "requires directory of bam files"}
OUT=${2?Error: "provide output directory name"}
for seed in 2 5 9
        do bash $SRC/batch_call_peak.sh ${DIR}/seed$seed ${OUT}
        done	
