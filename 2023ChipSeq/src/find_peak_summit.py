
import sys, os
import numpy as np
import pandas as pd

# input is .intersect.bed (from multiinter)
infile = sys.argv[1]
outfile = sys.argv[2]

def get_matching_peak(chrom, coordinate, file):
    # from one replicate (narrowPeak), list peaks within which a coordinate is contained
    # 18Apr2023 limit search in the same chrom
    df_replicate = pd.read_csv(file, sep='\t', header=None)
    matching_peak = []
    for index, row in df_replicate.iterrows():
        chromosome = row[0]
        start = row[1]
        end = row[2]
        offset = row[9]
        peak_height = row[6]
        summit = start + offset
        if chromosome == chrom:
            if start < coordinate < end:
                matching_peak.append([start, end, offset, peak_height])
    return matching_peak

with open(infile, 'r') as f: #get file name from df can cause error because df column names cannot be the same
    header = f.readline().strip().split('\t') 
replicates = header[5:] #changed from [-3:] to [5:] to allow only 2 replicate

df_intervals = pd.read_csv(infile, sep='\t')
common_intervals = []
for index, row in df_intervals.iterrows():
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    num = row['num'] # number of replicates where the interval is present
    mid = (start + end) // 2
    if num >= 2: # only proceed with common intervals present in >2 replicates
        offsets = []
        enrichment = []
        for rep in replicates:
            if row[rep] == 1: # if present
                matching_peak = get_matching_peak(chrom, mid, rep)
                if len(matching_peak) == 1: # only proceed when there is one match
                    offsets.append(matching_peak[0][2])
                    enrichment.append(matching_peak[0][3])
        avg_offset = int(np.mean(offsets))
        new_summit = start + avg_offset
        avg_enrichment = '{0:.2f}'.format(np.mean(enrichment))
        offsets = [str(i) for i in offsets]
        enrichment = ['{0:.2f}'.format(j) for j in enrichment]
        common_intervals.append([chrom, start, end, ';'.join(offsets), new_summit, ';'.join(enrichment), avg_enrichment])
header = ['chrom', 'start', 'end', 'original_peak_offset', 'average_summit', 'original_enrichment', 'average_enrichment']
df_common_intervals = pd.DataFrame(common_intervals, columns = header)
df_common_intervals.to_csv(outfile, sep='\t')
