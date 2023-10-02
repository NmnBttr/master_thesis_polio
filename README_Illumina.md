# Run CoVpipe2 for Polio Illumina Data

## SETUP
```bash
WORKDIR=/home/namuun/Documents/CoVpipe_Polio
FASTQ=/home/namuun/Documents/CoVpipe_Polio/data/230222_MS1_Run23-055/kreibichj/fastq
REF=/home/namuun/Documents/CoVpipe_Polio/Referenz_poliovirus_Sabin_1_converted.fasta
PRIMER=/home/namuun/Downloads/Poliovirus_1/results/polio/final_primers/0_aligned_consensus_masked/primer_polioS1.bed
# FASTQCONVERTED=/home/namuun/Documents/CoVpipe_Polio/data/230222_MS1_Run23-055/kreibichj/fastq/converted_data
```

## Login docker:
```bash
conda activate nextflow
sudo chmod 666 /var/run/docker.sock
```

# Prepare Reference
```bash
awk -F "|" '/^>/ {close(F); ID=$1; gsub("^>", "", ID); F=ID".fasta"} {print >> F}' Referenz_Poliovirus_Sabin_1_2_3.fasta

sed '/^[^>]/s/U/T/g' Referenz_poliovirus_Sabin_1.fasta > Referenz_poliovirus_Sabin_1_converted.fasta
REF=/home/namuun/Documents/CoVpipe_Polio/Referenz_poliovirus_Sabin_1_converted.fasta

# for f in Referenz_*; do 
#   sed '/^[^>]/s/u/t/g' | sed '/^[^>]/s/U/T/g' "$f" > "$f"; 
# done
```
# Run CoVpipe2
### Reshape primer bed file
```bash
awk '{print "AY184219", $2, $3, $4, $5, $6}' 0_aligned_consensus_masked.fasta_1000.primer.bed > primer_polio1.bed
awk '{if (NR%4==0 || NR%4==3){print $1"\t"$2"\t"$3"\t"$4"\tnCoV-2019_2\t"$6} else {print $1"\t"$2"\t"$3"\t"$4"\tnCoV-2019_1\t"$6}}' primer_polio1.bed > primer_polioS1.bed
```

### Convert FASTA REF file U to T
```bash
sed '/^[^>]/s/U/T/g' Referenz_poliovirus_Sabin_1.fasta > Referenz_poliovirus_Sabin_1_converted.fasta
for f in Referenz_*; do sed '/^[^>]/s/u/t/g' | sed '/^[^>]/s/U/T/g'; done
```
<!-- ### Convert FASTQ file U -> T 
```bash
gzip -cd 220223_23-01835_SiB_142_S1_L000_R1_001.fastq.gz | sed '/^[^>]/s/u/t/g' | sed '/^[^>]/s/U/T/g' | gzip > $FASTQCONVERTED/"test_file.fastq.gz"
for f in $FASTQ/*; do   
  cp "$f" "$f~" &&      
  gzip -cd "$f~" | sed '/^[^>]/s/u/t/g' | sed '/^[^>]/s/U/T/g' | gzip > "$f"; 
done

``` -->

### concat amplicon pools belonging to the same sample (Jonas Primer)

```bash
cat 220223_23-01840_SiB_147_S6_L000_R1_001.fastq.gz 220223_23-01841_SiB_148_S7_L000_R1_001.fastq.gz > SiB_147_148_S_6_7_L000_concat_R1_001.fastq.gz

cat 220223_23-01840_SiB_147_S6_L000_R2_001.fastq.gz 220223_23-01841_SiB_148_S7_L000_R2_001.fastq.gz > SiB_147_148_S_6_7_L000_concat_R2_001.fastq.gz

cat 220223_23-01842_SiB_149_S8_L000_R1_001.fastq.gz 220223_23-01843_SiB_150_S9_L000_R1_001.fastq.gz > SiB_149_150_S_8_9_L000_concat_R1_001.fastq.gz

cat 220223_23-01842_SiB_149_S8_L000_R2_001.fastq.gz 220223_23-01843_SiB_150_S9_L000_R2_001.fastq.gz > SiB_149_150_S_8_9_L000_concat_R2_001.fastq.gz

```

### Run nextflow for CoVpipe2
```bash
nextflow run rki-mf1/CoVpipe2 --ref_genome $REF --fastq /home/namuun/Documents/CoVpipe_Polio/test4/samplesheet.csv --list --cores 4 --max_cores 8 -r v0.4.2 -profile local,docker --primer_bed /home/namuun/Documents/CoVpipe_Polio/V_Primer_polio1/primer_polioS1.bed
# to continue 
# -resume
```

### Convert bam files to fastq to check if clipping worked
 
```bash
find *147*  -name "*.bam" > bam_to_fastq_jonas_primer.txt
find *149* -name "*.bam" >> bam_to_fastq_jonas_primer.txt 

cat bam_to_fastq_jonas_primer.txt | while read line; do samtools fastq $line > $line.fastq; done
```
### Get Matching reads to primer seqs

```bash

while read line ; do
  set $line
  # echo ${4^^} 
  grep ${4^^}$ SiB_151_152_23-0716_filtered.sorted.clipped.bam
  grep ^${4^^} SiB_151_152_23-0716_filtered.sorted.clipped.bam
done < <(tail -n +2 /home/namuun/Downloads/Poliovirus_1/results/polio/final_primers/0_aligned_consensus_masked/0_aligned_consensus_masked.fasta_1000.primer.tsv) > SiB_151_152_23-0716_filtered.sorted.clipped.bam.reads_with_primer.txt

while read line ; do
  set $line
  # echo ${4^^} 
  grep ${4^^}$ SiB_147_148_S_6_7.bam.fastq
  grep ^${4^^} SiB_147_148_S_6_7.bam.fastq
done < <(tail -n +2 /home/namuun/Downloads/Poliovirus_1/results/polio/final_primers/0_aligned_consensus_masked/0_aligned_consensus_masked.fasta_1000.primer.tsv) > SiB_147_148_S_6_7.bam.fastq.reads_with_primer.txt

# while read line ; do
#   set $line
#   seqkit grep --by-seq --max-mismatch 0 --pattern ${4^^} /home/namuun/Documents/CoVpipe_Polio/test6/results/02-Mapping/SiB_147_148_S_6_7/SiB_147_148_S_6_7.bam.fastq
# done < <(tail -n +2 /home/namuun/Downloads/Poliovirus_1/results/polio/final_primers/0_aligned_consensus_masked/0_aligned_consensus_masked.fasta_1000.primer.tsv) > SiB_147_148_S_6_7.bam.fastq.reads_with_primer.txt

# while read line ; do
#   set $line
#   seqkit grep --by-seq --max-mismatch 0 --pattern ${4^^} /home/namuun/Documents/CoVpipe_Polio/test6/results/02-Mapping/SiB_147_148_S_6_7/clipped/SiB_147_148_S_6_7.primerclipped.bam.fastq
# done < <(tail -n +2 /home/namuun/Downloads/Poliovirus_1/results/polio/final_primers/0_aligned_consensus_masked/0_aligned_consensus_masked.fasta_1000.primer.tsv) > SiB_147_148_S_6_7.primerclipped.bam.fastq.reads_with_primer.txt

# while read line ; do
#   set $line
#   seqkit grep --by-seq --max-mismatch 1 --pattern ${4^^} /home/namuun/Documents/CoVpipe_Polio/Illumina_complete_run_2/results/02-Mapping/SiB_149_150_S_8_9/SiB_149_150_S_8_9.bam.fastq
# done < <(tail -n +2 /home/namuun/Downloads/0_aligned_consensus_masked.fasta_1000.primer.tsv) > SiB_149_150_S_8_9.bam.fastq.reads_with_primer.txt

# while read line ; do
#   set $line
#   seqkit grep --by-seq --max-mismatch 1 --pattern ${4^^} /home/namuun/Documents/CoVpipe_Polio/Illumina_complete_run_2/results/02-Mapping/SiB_149_150_S_8_9/clipped/SiB_149_150_S_8_9.primerclipped.bam.fastq
# done < <(tail -n +2 /home/namuun/Downloads/0_aligned_consensus_masked.fasta_1000.primer.tsv) > SiB_149_150_S_8_9.primerclipped.bam.fastq.reads_with_primer.txt

```
<!-- PRIMER1="cgacyactgaattagccgcctc"
echo ${PRIMER1^^}
 seqkit grep --by-seq --max-mismatch 0 --pattern $4 /home/namuun/Documents/CoVpipe_Polio/Illumina_complete_run_2/results/02-Mapping/SiB_147_148_S_6_7/SiB_147_148_S_6_7.bam.fastq
seqkit grep --by-seq --max-mismatch 0 --pattern "CGCCTGTTTTATACTCCCTTCCC" -R 100:130 /home/namuun/Documents/CoVpipe_Polio/test6/results/02-Mapping/SiB_147_148_S_6_7/SiB_147_148_S_6_7.bam.fastq -->

### Get mutations from vcf file 

```bash
find . -name "*vcf.annotation.covered.af.vcf.gz" > vcf_files.txt
cat vcf_files.txt | while read line; do bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AB\n' $line  >> vcf_files.txt; done 

# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AB\n' SiB_146_S5.vcf.annotation.covered.af.vcf.gz
# AB = Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous
```

### Run FLASH
./flash $FASTQ/SiB_147_148_S_6_7_L000_concat_R1_001.fastq.gz $FASTQ/SiB_147_148_S_6_7_L000_concat_R2_001.fastq.gz -d $WORKDIR/FLASH_test -o "SiB_147_148" | tee $WORKDIR/FLASH_test/SiB_147_148_flash.log

./flash $FASTQ/SiB_149_150_S_8_9_L000_concat_R1_001.fastq.gz $FASTQ/SiB_149_150_S_8_9_L000_concat_R2_001.fastq.gz -d $WORKDIR/FLASH_test -o "SiB_149_150" | tee $WORKDIR/FLASH_test/SiB_149_150_flash.log

./flash $FASTQ/220223_23-01839_SiB_146_S5_L000_R1_001.fastq.gz $FASTQ/220223_23-01839_SiB_146_S5_L000_R2_001.fastq.gz -d $WORKDIR/FLASH_test -o "SiB_146" | tee $WORKDIR/FLASH_test/SiB_146_flash.log

./flash $FASTQ/220223_23-01838_SiB_145_S4_L000_R1_001.fastq.gz $FASTQ/220223_23-01838_SiB_145_S4_L000_R2_001.fastq.gz -d $WORKDIR/FLASH_test -o "SiB_145" | tee $WORKDIR/FLASH_test/SiB_145_flash.log

### sample Jonas primer samples so that 


#### See head of gzipped text file
zcat 220223_23-01835_SiB_142_S1_L000_R1_001.fastq.gz | head -n 5


nextflow run rki-mf1/CoVpipe2 --ref_genome $REF --fastq /home/namuun/Documents/CoVpipe_Polio/test4/samplesheet.csv --list --cores 4 --max_cores 8 -r v0.4.2 -profile local,docker --primer_bed $PRIMER


while read line ; do
  set $line
  echo ${4^^}
  #grep ${4^^}$ SiB_147_148_S_6_7.primerclipped.bam.fastq
done < <(tail -n +2 /home/namuun/Downloads/Poliovirus_1/results/polio/final_primers/0_aligned_consensus_masked/0_aligned_consensus_masked.fasta_1000.primer.tsv) 

grep ^AGCAATGGGAAAGAAGAAGAGAGA$ SiB_147_148_S_6_7.bam.fastq

CGCCTGTTTTATACTCCCTTCCC
CGACYACTGAATTAGCCGCCTC
ATGGGTGCTCAGGTTTCATCAC
GTCACATCAAATTCAGGCAGCG
TGCAATATTACCATTGGCCCCA
ATTTCCTTGGAGTGTGCTGGTC
TTTCGACACCCAGAGAGATGGA
ATCTTCCTGAGTGGCCAAATGG
TGTCCAGTTACGGAGGAAATTGG
TCRTAGGCATACAAGTCTCTAATGT
GTCACAGAATCAAGAGCCCAGG
TCAGGTCRTCCATAATCACCAC
AGAGCAAACACCGTATTGAACCA
TAGCCATAGCCACTGCRTAATC
TGTTYCAAGGACCACTCCAGTA
AGCACATTTGTTCTGTGTTGATGT
DGCCCTGAAGCGATCATACTTC
CATGGGGGTAGGAAGCAATTACA
AGCAATGGGAAAGAAGAAGAGAGA
TCCAATTCGACTGAGGTAGGGT
