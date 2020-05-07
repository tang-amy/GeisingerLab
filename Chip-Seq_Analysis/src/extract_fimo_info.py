

from sys import argv
import os, pandas as pd


DIR = argv[1] if len(argv)>1 else "/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/meme-chip"
df_out = pd.DataFrame(columns=['Filename', 'Start', 'End', 'IUPAC code', 'Match Sequence'])
for entry in os.listdir(DIR):
    full_path = os.path.join(DIR, entry)
    if os.path.isdir(full_path):
        path = os.path.join(full_path, 'bed_intersect_consensus_meme_chip', 'fimo_output')
        for filename in os.listdir(path):
            file = os.path.join(path, filename)
            df_fimo = pd.read_csv(file, sep='\t', skiprows=1,
                                  names=['chrom', 'fimo', 'motif', 'start', 'end', 'p_val', 'strand', '.', 'string'])
            summary = []
            for i in df_fimo.index.tolist():
                down_scale = file.split('/')[-4]
                name = os.path.splitext(os.path.basename(file))[0]
                start = df_fimo['start'][i]
                end = df_fimo['end'][i]
                strand = df_fimo['strand'][i]
                info = df_fimo['string'][i].split(';')
                IUPAC = info[0][5:]
                sequence = info[5][9:]
                summary.append([down_scale, name, start, end, strand, IUPAC, sequence])
            df_out = pd.DataFrame(summary,
                                  columns=['Downsample Scale', 'Filename', 'Start', 'End', 'Strand', 'IUPAC code', 'Match Sequence'])
            if not os.path.isfile('/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/meme-chip/fimo_summary.tsv'):
                df_out.to_csv('/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/meme-chip/fimo_summary.tsv',
                              mode='a', sep='\t', header='column_names', index=False)
            else:
                df_out.to_csv('/Volumes/Seagate/10022019_ChipSeq/mapped_SAM/meme-chip/fimo_summary.tsv',
                              mode='a', sep='\t', header=False, index=False)
