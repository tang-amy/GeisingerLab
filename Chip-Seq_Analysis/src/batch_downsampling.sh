## Yunfei Dai
## 04/28/2020

"""
downsample all .bam files in directory ./sorted_bam/
for each file, 3 downsampled bam files containing ~100,000 reads are created, using 3 seed numbers (2,5,9) 
exact read number from each downsampled files are counted and writted into downsampling_log.txt 
"""

for file in sorted_bam/*.bam
do fname=$(basename $file); 
read_count=$(samtools view $file | wc -l); 
scale=$(echo "scale=6; 100000/$read_count" | bc);
scale2=$(echo "scale=6; 2+$scale" | bc);
scale5=$(echo "scale=6; 5+$scale" | bc);
scale9=$(echo "scale=6; 9+$scale" | bc);
ofile2=downsampled_bam/seed2/${fname/sorted.bam/s2_100k.sorted.bam};
ofile5=downsampled_bam/seed5/${fname/sorted.bam/s5_100k.sorted.bam};
ofile9=downsampled_bam/seed9/${fname/sorted.bam/s9_100k.sorted.bam};
samtools view -s $scale2 -b $file -o $ofile2;
samtools view -s $scale5 -b $file -o $ofile5;
samtools view -s $scale9 -b $file -o $ofile9;
echo $fname >> downsampling_log.txt;
echo "s2_100k: $(samtools view $ofile2 | wc -l) reads" >> downsampling_log.txt;
echo "s5_100k: $(samtools view $ofile5 | wc -l) reads" >> downsampling_log.txt;
echo "s9_100k: $(samtools view $ofile9 | wc -l) reads" >> downsampling_log.txt;
done
