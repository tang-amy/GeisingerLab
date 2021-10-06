## Yunfei Dai
## 10/06/2021

from sys import argv
from Bio import SeqIO, Entrez
from matplotlib import pyplot as plt

try:
    infile = argv[1]
except Exception:
    print("Please provide input table.")

try:
    outfile = argv[2]
except Exception:
    outfile = "output_histogram.pdf"

    
plt.figure(figsize=(8,6))
plt.hist(all_seq_lengths, bins='auto', alpha=0.5, label="all sequences")
plt.hist(CDD_hit_length, bins='auto', alpha=0.5, label="subset_CDD_and_SigIP")
plt.hist(predisi_TMHMM_hit_length, bins='auto', alpha=0.5, label="subset_CDD_SigIP_TMHMM_and_predisi")
plt.hist(size_excluded_lengths, bins='auto', alpha=0.5, label="size_excluded (final shortlist)")

plt.xlabel("Sequence Length")
plt.ylabel("Count")
plt.legend(loc="upper right")
plt.show()
