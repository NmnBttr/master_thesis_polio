
```bash
REF=/home/namuun/Documents/Polio_Nanopore_Analysis/Referenz_poliovirus_Sabin_1_converted.fasta
BAM=/home/namuun/Downloads/20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz.sorted.bam.clipped.sort.bam
FASTQ=/home/namuun/Downloads/20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.fastq.gz.filtered.fastq

#1
samtools dict $REF -o $(basename $REF).dict
#2
java -jar /home/namuun/jvarkit/dist/jvarkit.jar sortsamrefname $BAM -o $(basename $BAM).refsorted.bam
#3
java -jar /home/namuun/jvarkit/dist/jvarkit.jar biostar154220 -n 500 --samoutputformat BAM -o $(basename $BAM).coverage.bam --regions '[ 743\t1856\n;,2370\t3463\n;,4467\t5472\n;,5953\t6996\n;]'
20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz.sorted.bam.clipped.sort.bam.refsorted.bam
#4 
samtools sort -O BAM --reference $REF 20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz.sorted.bam.clipped.sort.bam.coverage.bam > 20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz.sorted.bam.clipped.sort.bam.coverage.sorted.bam
#5
samtools index -b 20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz.sorted.bam.clipped.sort.bam.coverage.sorted.bam
#6 
samtools fastq 20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz.sorted.bam.clipped.sort.bam.coverage.sorted.bam > 20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz.sorted.bam.clipped.sort.bam.coverage.sorted.bam.fastq
#7 - FLYE
flye --meta -g 7400 --nano-corr 20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz.sorted.bam.clipped.sort.bam.coverage.sorted.bam.fastq -o ./ --min-overlap 1000
#7.2 - FLYE ohne bam capping
flye --meta -g 7400 --nano-raw /home/namuun/Downloads/20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz -o ./ --min-overlap 1000
#8 plot coverage distribution over the genome
bamdash -b 20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz.sorted.bam.clipped.sort.bam.coverage.sorted.bam -r AY184219
#9 de novo assembly with haploflow
haploflow --read-file 20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz.sorted.bam.clipped.sort.bam.coverage.sorted.bam.fastq --out haploflow-cap-500 --debug 1

#10 de novo assembly illumina merge R1 and R2 of illumina data
conda install bioconda::fastq-join
FASTQ_R1=/home/namuun/Downloads/230418_AS-2882_SiB179_S1_L000_R1_001.fastq.gz
FASTQ_R2=/home/namuun/Downloads/230418_AS-2882_SiB179_S1_L000_R2_001.fastq.gz
fastq-join $FASTQ_R1 $FASTQ_R2 -o 230418_AS-2882_SiB179_S1_L000_.fastq

conda activate haploflow
haploflow --read-file 230418_AS-2882_SiB179_S1_L000_.fastq --out $PWD --debug 1

conda activate bioawk
bioawk -c fastx '{print $name, length($seq)}' < contigs.fa 


```
