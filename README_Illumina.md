# Run CoVpipe2 for Polio Illumina Data

## SETUP
```bash
srun --cpus-per-task=4 --mem=50GB --gres=local:30 --gpus=1 --pty --time=03:00:00 bash -i 

WORKDIR=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/CoVpipe2_run_updated_primer
FASTQ=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/data/2023-08-18-PolioPrimer-Illumina-runID115
REF=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/References/Referenz_Poliovirus_Sabin_1.fasta
REF2=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/References/Referenz_Poliovirus_Sabin_2.fasta
PRIMER=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/primer/polio1/primer_polioS1.bed
PRIMER2=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/primer/polio2/primer_polioS2.bed
REF3=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/References/Referenz_Poliovirus_Sabin_3.fasta
PRIMER3=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/primer/polio3/primer_polioS3.bed
```

## Login docker (local only):
```bash
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

### Update primer positions in bed file according to the Reference
```bash
conda env create -f py_env_bed_pos_correction.yml
conda activate py_bed
# tsv and bed files are the direct output from varvamp
python correct_primer_positions.py $REF $PWD/polio1/*.primer.tsv $PWD/polio1/*.primer.bed
```
### Reshape primer bed file
```bash
awk '{print "AY184219""\t", $2"\t", $3"\t", $4"\t", $5"\t", $6"\t"}' *position.corrected.bed > primer_polio1.bed
awk '{if (NR%4==0 || NR%4==3){print $1"\t"$2"\t"$3"\t"$4"\tnCoV-2019_2\t"$6} else {print $1"\t"$2"\t"$3"\t"$4"\tnCoV-2019_1\t"$6}}' primer_polio1.bed > primer_polioS1.bed

awk '{print "AY184220""\t", $2"\t", $3"\t", $4"\t", $5"\t", $6"\t"}' *position.corrected.bed > primer_polio2.bed
awk '{if (NR%4==0 || NR%4==3){print $1"\t"$2"\t"$3"\t"$4"\tnCoV-2019_2\t"$6} else {print $1"\t"$2"\t"$3"\t"$4"\tnCoV-2019_1\t"$6}}' primer_polio2.bed > primer_polioS2.bed

awk '{print "AY184221""\t", $2"\t", $3"\t", $4"\t", $5"\t", $6"\t"}' *position.corrected.bed > primer_polio3.bed
awk '{if (NR%4==0 || NR%4==3){print $1"\t"$2"\t"$3"\t"$4"\tnCoV-2019_2\t"$6} else {print $1"\t"$2"\t"$3"\t"$4"\tnCoV-2019_1\t"$6}}' primer_polio3.bed > primer_polioS3.bed
```

### Convert FASTA REF file U to T
```bash
sed '/^[^>]/s/U/T/g' Referenz_poliovirus_Sabin_1.fasta > Referenz_poliovirus_Sabin_1_converted.fasta
for f in Referenz_*; do sed '/^[^>]/s/u/t/g' | sed '/^[^>]/s/U/T/g'; done
```

### concat amplicon pools belonging to the same sample (Jonas Primer)

```bash
cat 220223_23-01840_SiB_147_S6_L000_R1_001.fastq.gz 220223_23-01841_SiB_148_S7_L000_R1_001.fastq.gz > SiB_147_148_S_6_7_L000_concat_R1_001.fastq.gz

cat 220223_23-01840_SiB_147_S6_L000_R2_001.fastq.gz 220223_23-01841_SiB_148_S7_L000_R2_001.fastq.gz > SiB_147_148_S_6_7_L000_concat_R2_001.fastq.gz

cat 220223_23-01842_SiB_149_S8_L000_R1_001.fastq.gz 220223_23-01843_SiB_150_S9_L000_R1_001.fastq.gz > SiB_149_150_S_8_9_L000_concat_R1_001.fastq.gz

cat 220223_23-01842_SiB_149_S8_L000_R2_001.fastq.gz 220223_23-01843_SiB_150_S9_L000_R2_001.fastq.gz > SiB_149_150_S_8_9_L000_concat_R2_001.fastq.gz

```
## Install nextflow for CovPipe2:
```bash
conda create -n nextflow -c bioconda nextflow
conda activate nextflow
```

### Run nextflow for CoVpipe2
```bash
cd run_Polio1_2
nextflow run rki-mf1/CoVpipe2 --ref_genome $REF --fastq /scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/CoVpipe2_run_updated_primer/Polio1/samplesheet.csv --list -r v0.4.3 -profile slurm,singularity --primer_bed $PRIMER --singularity_cache_dir /biolibs/mf1/singularity/CoVpipe2/ 

cd run_Polio2
nextflow run rki-mf1/CoVpipe2 --ref_genome $REF2 --fastq $PWD/samplesheet.csv --list -r v0.4.3 -profile slurm,singularity --primer_bed $PRIMER2 --singularity_cache_dir /biolibs/mf1/singularity/CoVpipe2/ 

cd run_Polio3
nextflow run rki-mf1/CoVpipe2 --ref_genome $REF3 --fastq $PWD/samplesheet.csv --list -r v0.4.3 -profile slurm,singularity --primer_bed $PRIMER3 --singularity_cache_dir /biolibs/mf1/singularity/CoVpipe2/

# to continue 
# -resume
# nextflow run rki-mf1/CoVpipe2 -r v0.4.3 -profile slurm,singularity,test --singularity_cache_dir /biolibs/mf1/singularity/CoVpipe2/

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


### Get mutations from vcf file 

```bash
find ./Polio*/results/03-Variant-Calling/ -name "*vcf.annotation.covered.af.vcf.gz" > vcf_files.txt
# cat vcf_files.txt | while read line; do bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AB\n' $line  >> vcf_files.txt; done 
while read line; do
  # bcftools view --min-af 0.1 --exclude-types indels $line > $WORKDIR/filtered_vcf/$(basename $line ).filtered.vcf
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%DP\t[%AD]\n' $line  > $PWD/vcf_snps/$(basename $line).txt
done < vcf_files.txt

# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AB\n' SiB_146_S5.vcf.annotation.covered.af.vcf.gz
# AB = Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous
```

### Run FLASH -> https://github.com/genome-vendor/FLASH/blob/master/MANUAL
```bash
./flash $FASTQ/SiB_147_148_S_6_7_L000_concat_R1_001.fastq.gz $FASTQ/SiB_147_148_S_6_7_L000_concat_R2_001.fastq.gz -d $WORKDIR/FLASH_test -o "SiB_147_148" | tee $WORKDIR/FLASH_test/SiB_147_148_flash.log

./flash $FASTQ/SiB_149_150_S_8_9_L000_concat_R1_001.fastq.gz $FASTQ/SiB_149_150_S_8_9_L000_concat_R2_001.fastq.gz -d $WORKDIR/FLASH_test -o "SiB_149_150" | tee $WORKDIR/FLASH_test/SiB_149_150_flash.log

./flash $FASTQ/220223_23-01839_SiB_146_S5_L000_R1_001.fastq.gz $FASTQ/220223_23-01839_SiB_146_S5_L000_R2_001.fastq.gz -d $WORKDIR/FLASH_test -o "SiB_146" | tee $WORKDIR/FLASH_test/SiB_146_flash.log

./flash $FASTQ/220223_23-01838_SiB_145_S4_L000_R1_001.fastq.gz $FASTQ/220223_23-01838_SiB_145_S4_L000_R2_001.fastq.gz -d $WORKDIR/FLASH_test -o "SiB_145" | tee $WORKDIR/FLASH_test/SiB_145_flash.log
```

### Run BAMdash for Coverage Plots 
```bash
mkdir run_polio1
cd run_polio1
find $WORKDIR/Polio1/results -name "*.primerclipped.bam" > bam_files.txt
while read line; do
  bamdash -b $line -r AY184219 && mv *_plot.html $(basename $line).html
done < bam_files.txt

mkdir run_polio2
cd run_polio2
find $WORKDIR/Polio2/results -name "*.primerclipped.bam" > bam_files.txt
while read line; do
  bamdash -b $line -r AY184220 && mv *_plot.html $(basename $line).html
done < bam_files.txt

mkdir run_polio3
cd run_polio3
find $WORKDIR/Polio3/results -name "*.primerclipped.bam" > bam_files.txt
while read line; do
  bamdash -b $line -r AY184221 && mv *_plot.html $(basename $line).html
done < bam_files.txt

# bamdash -b SiB179.bam -r AY184219
```
#### See head of gzipped text file
```bash
zcat 220223_23-01835_SiB_142_S1_L000_R1_001.fastq.gz | head -n 5
```

# de novo assembly using SPAdes
```bash
conda activate spades 

spades.py --rnaviral -o SiB194 -1 $FASTQ/230418_AS-2897_SiB194_S16_L000_R1_001.fastq.gz -2 $FASTQ/230418_AS-2897_SiB194_S16_L000_R2_001.fastq.gz


spades.py --rnaviral -o SiB195 -1 $FASTQ/230418_AS-2898_SiB195_S17_L000_R1_001.fastq.gz -2 $FASTQ/230418_AS-2898_SiB195_S17_L000_R2_001.fastq.gz
spades.py --rnaviral -o SiB196 -1 $FASTQ/230418_AS-2899_SiB196_S18_L000_R1_001.fastq.gz -2 $FASTQ/230418_AS-2899_SiB196_S18_L000_R2_001.fastq.gz
spades.py --rnaviral -o SiB197 -1 $FASTQ/230418_AS-2900_SiB197_S19_L000_R1_001.fastq.gz -2 $FASTQ/230418_AS-2900_SiB197_S19_L000_R2_001.fastq.gz
spades.py --rnaviral -o SiB198 -1 $FASTQ/230418_AS-2901_SiB198_S20_L000_R1_001.fastq.gz -2 $FASTQ/230418_AS-2901_SiB198_S20_L000_R2_001.fastq.gz
spades.py --rnaviral -o SiB199 -1 $FASTQ/230418_AS-2902_SiB199_S21_L000_R1_001.fastq.gz -2 $FASTQ/230418_AS-2902_SiB199_S21_L000_R2_001.fastq.gz

# Rekombinanten
#188 2891_SiB188_S10
#194 2897_SiB194_S16
#195 2898_SiB195_S17
#196 2899_SiB196_S18
#197 2900_SiB197_S19
#198 2901_SiB198_S20
#199 2902_SiB199_S21

conda activate bioawk
for READS in $PWD/SiB*/contigs.fasta; do
    echo $READS 
    bioawk -c fastx '{ print $name, length($seq) }' | sort -k1,1rn | head -1 < $READS 
done > longest_contig_sizes.txt


for READS in $PWD/SiB*/contigs.fasta; do
  echo $READS 
  cat $READS | bioawk -c fastx '{ print length($seq), $name }' | sort -k1,1rn | head -1 
done > longest_contig_sizes.tx

```
## get consensus files
```bash
for READS in $PWD/SiB*/*.iupac_consensus.fasta; do
  cp $READS /scratch/Tausch/batturn/consensus_sequences_polio3
done
```