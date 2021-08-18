
import csv
CDD_list = "/Users/yunfei/Downloads/GT9TM664013_CDDhits_standard.txt"
outfile = "/Users/yunfei/Downloads/GT9TM664013_CDDhits_standard_from_py.csv"
with open(CDD_list, 'r') as CDD:
    writer_list = []
    for line in CDD:
        line = line.split("\t")
        line = [i.strip('\"') for i in line]
        line = [i.strip() for i in line]
        print(line)
        if len(line) > 1 and len(str(line)) > 1:
            writer_list.append(line)

with open(outfile, 'w') as output:
    writer = csv.writer(output, delimiter='\t')
    writer.writerows(writer_list)
