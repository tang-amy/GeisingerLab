
from sys import argv
import os, codecs , pandas as pd


file = "/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/meme-chip/1M_downsampled/" \
       "bed_intersect_consensus_meme_chip/html_output/BfmR-ChIP-28_seed9.consensus_peak.meme-chip.html"

f = codecs.open(file, 'r', 'utf-8')
linelist = list(enumerate(f, 1))

for num, line in linelist:
       path = os.path.abspath(file)
       size = path.split('/')[-4][:-12]
       result_name = path.split('/')[-2][:-15] + '_consensus'
       seed = path.split('/')[-2][-16:-15]
       INFO = []
       if "\"consensus\":" in line:
              method = linelist[num-2][1].strip().split("\"")[-2]
              motif = linelist[num-1][1].strip().split("\"")[-2]
              nsites = linelist[num+2][1].strip().split(": ")[-1][:-1]
              e_value = float(linelist[num+3][1].strip().split("\"")[-2])
              print(method)
              print(motif)
              print(nsites)
              print(e_value)

