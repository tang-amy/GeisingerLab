from Bio import SeqIO

infile = '/Users/yunfei/Downloads/JBA71FLAG.consensus_peak.fasta'
reference = '/Users/yunfei/GeisingerLab/CRISPRi/reference_files/NZ_CP012004.gbk'

with open(infile) as f:
    for line in f:
        if '>' in line:
            line = line.strip().split(':')[1]
            coordinate = line.split('-')
            coord_up = coordinate[0]
            coord_down = coordinate[1]
            #print(coord_up+','+coord_down)

def gene_dic(gbk):
    # Parse gbk file to get dictionaries containing gene positions and locus tags
    recs = [rec for rec in SeqIO.parse(gbk, "genbank")]
    gbk_dic = {}
    for rec in recs:
        feats = [feat for feat in rec.features if feat.type == "CDS"]
        for feat in feats:
            if 'protein_id' in feat.qualifiers:
                protein = ''.join(feat.qualifiers['protein_id'])
                tag = ''.join(feat.qualifiers['locus_tag'])
                gbk_dic.update({protein: tag})
    return gbk_dic