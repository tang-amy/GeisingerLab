


from optparse import OptionParser
import pandas as pd


options = OptionParser()
options.add_option("-i", "--tss_input", dest="tss",
                   help="input list of TSS sites)")
options.add_option("-o", "--tss_output", dest="outfile",
                   help="input list of TSS sites)")


def primary_tss(tss, outfile):
    df_TSS = pd.read_csv(tss, sep='\t')
    tss_index = df_TSS.index.tolist()[:-1]
    tss_coordinate = [int(i) for i in df_TSS["TSS coordinate"].tolist()[:-1]]
    tss_type = [int(i) for i in df_TSS["Primary"].tolist()[:-1]]
    tss_strand = [i for i in df_TSS["Strand"].tolist()[:-1]]
    tss_tag = [i for i in df_TSS["Locus_tag"].tolist()[:-1]]
    primary_tss = []
    counter = 0
    for index in tss_index:
        if tss_type[index] != 0:
            primary_tss.append([tss_coordinate[index], tss_tag[index], tss_strand[index]])
            counter += 1

    df = pd.DataFrame(primary_tss,
                      columns=['TSS coordinate', 'Locus_tag', 'Strand'])
    df.to_csv(outfile, sep='\t', index=False)
    print('Found ', counter, ' primary TSSs out of ', tss_index[-1], ' records.')


def main():
    opts, args = options.parse_args()
    tss = opts.tss
    outfile = opts.outfile
    primary_tss(tss, outfile)


if __name__ == '__main__':
    main()
