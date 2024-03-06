# SETUP

```bash
srun --cpus-per-task=4 --mem=50GB --gres=local:30 --gpus=1 --pty --time=03:00:00 bash -i

WORKDIR=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/Nanopore_analysis_updated_primer
FASTQ=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/data/2023-08-18-PolioPrimer-ONT-runID114

REF=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/References/Referenz_Poliovirus_Sabin_1.fasta
REF2=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/References/Referenz_Poliovirus_Sabin_2.fasta
REF3=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/References/Referenz_Poliovirus_Sabin_3.fasta
REF123=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/References/Referenz_Poliovirus_Sabin_1_2_3.fasta
PRIMER=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/primer/polio1/primer_polioS1.bed
PRIMER2=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/primer/polio2/primer_polioS2.bed
PRIMER3=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/primer/polio3/primer_polioS3.bed

```

# unzip fastq files
```bash
#!/bin/sh
for zip in *fastq_pass.zip
do
  dirname=`echo $zip | sed 's/\.zip$/_barcode0/'`
  if mkdir "$dirname"
  then
    if cd "$dirname"
    then
      unzip ../"$zip"
      cd ..
      # rm -f $zip # Uncomment to delete the original zip file
    else
      echo "Could not unpack $zip - cd failed"
    fi
  else
    echo "Could not unpack $zip - mkdir failed"
  fi
done
```

# prepare reads,  samples
```bash
for READS in $FASTQ/*fastq_pass; do BN=$(basename $READS); cat $READS/*.fastq.gz > $BN.fastq.gz; done
```


# prepare reads, 2 samples each 2 seqs for preliminary analysis
```bash
for READS in $FASTQ/barcode*; do BN=$(basename $READS); cat $READS/*.fastq.gz > $BN.fastq.gz; done

cat barcode21_SiB_151_23-01715.fastq.gz barcode22_SiB_152_23-01716.fastq.gz > SiB_151_152_23-0716.fastq.gz
cat barcode23_SiB_153_23-01717.fastq.gz barcode24_SiB_154_23-01718.fastq.gz > SiB_153_154_23-01717.fastq.gz
```

# Filter reads for the expected amplicon length (+/- 100 bp)

```bash
cd Filtlong
conda activate filtlong

for READS in $WORKDIR/fastq_barcode_merged/*.fastq.gz; do
    bin/filtlong -a $REF --min_length 100 --keep_percent 90 --trim $READS | gzip > $WORKDIR/fastq_barcode_merged/$(basename $READS .fastq.gz).filtered.fastq.gz
done

#  Anylysis run with corrected bed positions
RAWFASTQ=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/Nanopore_analysis/fastq_barcode_merged

for READS in $RAWFASTQ/run_polio1/*fastq_pass.fastq.gz; do
    bin/filtlong -a $REF --min_length 500 --keep_percent 90 --trim $READS | gzip > $WORKDIR/data_preprocessing/Polio1/$(basename $READS).filtered.fastq.gz
done

for READS in $RAWFASTQ/run_polio2/*fastq_pass.fastq.gz; do
    bin/filtlong -a $REF2 --min_length 500 --keep_percent 90 --trim $READS | gzip > $WORKDIR/data_preprocessing/Polio2/$(basename $READS).filtered.fastq.gz
done

for READS in $RAWFASTQ/run_polio3/*fastq_pass.fastq.gz; do
    bin/filtlong -a $REF3 --min_length 500 --keep_percent 90 --trim $READS | gzip > $WORKDIR/data_preprocessing/Polio3/$(basename $READS).filtered.fastq.gz
done

### Filtering with Chopper to try out filtering with different QC scores in order to check if it impacts any of the known coverage drops such as ~3k region and ~5.5-6k region in PV1
conda activate chopper

for READS in $RAWFASTQ/run_polio1/*fastq_pass.fastq.gz; do
    gunzip -c $READS | chopper -q 20 -l 500 | gzip > $WORKDIR/data/Polio1/$(basename $READS).filtered.fastq.gz
done

for READS in $RAWFASTQ/run_polio2/*fastq_pass.fastq.gz; do
    gunzip -c $READS | chopper -q 20 -l 500 | gzip > $WORKDIR/data/Polio2/$(basename $READS).filtered.fastq.gz
done

for READS in $RAWFASTQ/run_polio3/*fastq_pass.fastq.gz; do
    gunzip -c $READS | chopper -q 20 -l 500 | gzip > $WORKDIR/data/Polio3/$(basename $READS).filtered.fastq.gz
done
```


# Map the filtered reads to a reference genome
# before this step manually sort samples into separate directories for different references
```bash
conda activate minimap2
cd /scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/Nanopore_analysis/minimap2/minimap2-2.26_x64-linux

# Analysis run with updated Primer 
for READS in $WORKDIR/data_preprocessing/Polio1/*.filtered.fastq.gz; do
    ./minimap2 -ax map-ont $REF $READS > $WORKDIR/data_preprocessing/Polio1/$(basename $READS).sam
    samtools view -bS $WORKDIR/data_preprocessing/Polio1/$(basename $READS).sam | samtools sort - -@ 2 -o $WORKDIR/data_preprocessing/Polio1/$(basename $READS).sorted.bam
    samtools index $WORKDIR/data_preprocessing/Polio1/$(basename $READS).sorted.bam
done

for READS in $WORKDIR/data_preprocessing/Polio2/*.filtered.fastq.gz; do
    ./minimap2 -ax map-ont $REF2 $READS > $WORKDIR/data_preprocessing/Polio2/$(basename $READS).sam
    samtools view -bS $WORKDIR/data_preprocessing/Polio2/$(basename $READS).sam | samtools sort - -@ 2 -o $WORKDIR/data_preprocessing/Polio2/$(basename $READS).sorted.bam
    samtools index $WORKDIR/data_preprocessing/Polio2/$(basename $READS).sorted.bam
done

for READS in $WORKDIR/data_preprocessing/Polio3/*.filtered.fastq.gz; do
    ./minimap2 -ax map-ont $REF3 $READS > $WORKDIR/data_preprocessing/Polio3/$(basename $READS).sam
    samtools view -bS $WORKDIR/data_preprocessing/Polio3/$(basename $READS).sam | samtools sort - -@ 2 -o $WORKDIR/data_preprocessing/Polio3/$(basename $READS).sorted.bam
    samtools index $WORKDIR/data_preprocessing/Polio3/$(basename $READS).sorted.bam
done

```

# Primer clippping for preliminary analysis

```bash
conda activate sambcfenv

for BAMF in $PWD/*.sorted.bam; do
    samtools ampliconclip -b $PRIMER2 --strand --both-ends $BAMF -o $(basename $BAMF).clipped.bam | tee $(basename $BAMF).samtools.log
    samtools sort $(basename $BAMF).clipped.bam > $(basename $BAMF).clipped.sort.bam
    samtools index $(basename $BAMF).clipped.sort.bam
    samtools fastq $(basename $BAMF).clipped.sort.bam > $(basename $BAMF).clipped.sort.bam.fastq
done

for BAMF in $PWD/*.sorted.bam; do
    samtools ampliconclip -b $PRIMER2 --strand --both-ends --clipped $BAMF -o $PWD/ampliconclip/$(basename $BAMF).clipped.bam | tee $PWD/ampliconclip/$(basename $BAMF).samtools.log
    samtools sort $PWD/ampliconclip/$(basename $BAMF).clipped.bam > $PWD/ampliconclip/$(basename $BAMF).clipped.sort.bam
    samtools index $PWD/ampliconclip/$(basename $BAMF).clipped.sort.bam
    samtools fastq $PWD/ampliconclip/$(basename $BAMF).clipped.sort.bam > $PWD/ampliconclip/$(basename $BAMF).clipped.sort.bam.fastq
done 
# bamclipper.sh -b /home/namuun/Documents/Polio_Nanopore_Analysis/SiB_151_152_23-0716_filtered.sorted.bam -p $PRIMER
```

# Primer clipping with iVAR for 2nd analysis with updated primer 
```bash
conda activate ivar

for BAMF in $PWD/*.sorted.bam; do     
    ivar trim -b $PRIMER -p $(basename $BAMF).clipped -m 0 -q 0 -i $BAMF 
    samtools sort $(basename $BAMF).clipped.bam > $(basename $BAMF).clipped.sort.bam
    samtools index $(basename $BAMF).clipped.sort.bam
done
cd Polio2
for BAMF in $PWD/*.sorted.bam; do     
    ivar trim -b $PRIMER2 -p $(basename $BAMF).clipped -m 0 -q 0 -i $BAMF 
    samtools sort $(basename $BAMF).clipped.bam > $(basename $BAMF).clipped.sort.bam
    samtools index $(basename $BAMF).clipped.sort.bam
done
cd ..
cd Polio3
for BAMF in $PWD/*.sorted.bam; do     
    ivar trim -b $PRIMER3 -p $(basename $BAMF).clipped -m 0 -q 0 -i $BAMF 
    samtools sort $(basename $BAMF).clipped.bam > $(basename $BAMF).clipped.sort.bam
    samtools index $(basename $BAMF).clipped.sort.bam
done

```

#### check primer clipping
```bash
for BAMF in $PWD/*.fastq; do
    while read line ; do
        set $line 
        grep ${4^^}$ $BAMF
        grep ^${4^^} $BAMF
    done < <(tail -n +2 /home/namuun/Downloads/Poliovirus_1/results/polio/final_primers/0_aligned_consensus_masked/0_aligned_consensus_masked.fasta_1000.primer.tsv) > $(basename $BAMF).txt
done
```

# Run BAMdash for Coverage Plots 
## -r AY184219 for Polio1; -r AY184220 for Polio 2; -r AY184221 for Polio 3
```bash
conda activate BAMdash

# Analysis 2nd run with updated primer
for BAMF in $WORKDIR/data_preprocessing/Polio1/*.clipped.sort.bam; do
  bamdash -b $BAMF -r AY184219 && mv *_plot.html $(basename $BAMF).html 
done
cd ..
cd Polio2
for BAMF in $WORKDIR/data_preprocessing/Polio2/*.clipped.sort.bam; do
  bamdash -b $BAMF -r AY184220 && mv *_plot.html $(basename $BAMF).html 
done
cd ..
cd Polio3
for BAMF in $WORKDIR/data_preprocessing/Polio3/*.clipped.sort.bam; do
  bamdash -b $BAMF -r AY184221 && mv *_plot.html $(basename $BAMF).html 
done

# Piranha BAM
PIRANHA_OUTDIR=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/Nanopore_analysis/piranha_run/analysis_with_ref123_2023-11-09
for line in $PIRANHA_OUTDIR/barcode*/reference_analysis/AY184221/medaka_haploid_variant/*.bam; do
    html_name=$(echo  $(dirname $line) | cut -d'/' -f9)
    tab3=$(echo  $(dirname $line) | cut -d'/' -f11)
    bamdash -b $line -r AY184221 && mv *_plot.html "$html_name"_"$tab3".html"" 
done
```

# Variant calling with basecall_model_version_id=dna_r1041_e82_260bps_sup@v3.5.2
# Run Medaka
# ATTENTION: it is always good to assign an appropriate Medaka model based on the performed basecalling! 
```bash
conda activate medaka
for READS in $WORKDIR/fastq_barcode_merged/run_polio3/*.sort.bam; do
    medaka consensus $READS --model r1041_e82_400bps_sup_v4.0.0 --threads 2 $(basename $READS .output_consensus).hdf;
    medaka snp $REF3 $(basename $READS .output_consensus).hdf $(basename $READS ).snp.vcf --verbose --threshold 0.1;       
done

for READS in $PWD/*.sort.bam; do
    medaka consensus $READS --model r1041_e82_400bps_sup_v4.0.0 --threads 2 $(basename $READS .output_consensus).hdf;
    medaka snp $REF2 $(basename $READS .output_consensus).hdf $(basename $READS ).snp.vcf --verbose --threshold 0.1;       
done

### 2nd run with corrected primer bed 
for READS in $WORKDIR/data_preprocessing/Polio2/*.sort.bam; do
    medaka_consensus -i $READS -d $REF2 -m r1041_e82_400bps_sup_v4.0.0 -t 2 -o $(basename $READS);
    # medaka snp $REF $(basename $READS).hdf $(basename $READS).snp.vcf --verbose --threshold 0.1;       
done

for READS in $WORKDIR/data_preprocessing/Polio2/*.sort.bam; do
    medaka consensus $READS --model r1041_e82_400bps_sup_v4.0.0 --threads 2 $(basename $READS).hdf;
    medaka snp $REF2 $(basename $READS).hdf $(basename $READS).snp.vcf --verbose --threshold 0.1;       
done

for READS in $WORKDIR/data_preprocessing/Polio3/*.sort.bam; do
    medaka consensus $READS --model r1041_e82_400bps_sup_v4.0.0 --threads 2 $(basename $READS).hdf;
    medaka snp $REF3 $(basename $READS).hdf $(basename $READS).snp.vcf --verbose --threshold 0.1;       
done

#### medaka haploid variant

for READS in $WORKDIR/fastq_barcode_merged/run_polio1/*.sort.bam; do
    medaka_haploid_variant -i $READS -r $REF -o $(basename $READS) -f -x
done
### 2nd run with corrected bed file
for READS in $WORKDIR/data_preprocessing/Polio1/*.sort.bam; do
    medaka_haploid_variant -i $READS -r $REF -o $(basename $READS) -f -x
done
cd ..
cd Polio2
for READS in $WORKDIR/data_preprocessing/Polio2/*.sort.bam; do
    medaka_haploid_variant -i $READS -r $REF2 -o $(basename $READS) -f -x
done
cd ..
cd Polio3
for READS in $WORKDIR/data_preprocessing/Polio3/*.sort.bam; do
    medaka_haploid_variant -i $READS -r $REF3 -o $(basename $READS) -f -x
done

#### medaka variant -> already deprecated, better use clair3 
```
## Freebayes variant calling

```bash
for READS in $WORKDIR/samtools_ampliconclip_both_ends/*.sort.bam; do
    freebayes -f $REF $READS --min-alternate-count 10 --min-alternate-fraction 0.1 --min-coverage 20 --pooled-continuous --no-indels --no-mnps --haplotype-length -1  -v $(basename $READS).vcf;
done
for READS in $PWD/*.vcf; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AB\n' $READS > $(basename $READS ).txt
done
```

## NanoCaller variant calling

```bash
#sudo chmod 666 /var/run/docker.sock

#docker run --rm -it -v $PWD:$PWD -w $PWD genomicslab/nanocaller:3.4.1 /bin/bash

for READS in $WORKDIR/fastq_barcode_merged/run_polio1/*.sort.bam; do
    NanoCaller --bam $READS --ref $REF --preset ont --mode snps --mincov 20 --min_allele_freq 0.1 --output $PWD/run_polio1 --prefix $(basename $READS)
done

for READS in $PWD/*.sort.bam; do
    NanoCaller --bam $READS --ref $REF2 --preset ont --mode snps --mincov 20 --min_allele_freq 0.1 --output $PWD --prefix $(basename $READS).nanocaller
done

### 2nd run with corrected bed
for READS in $WORKDIR/data_preprocessing/Polio1/*.sort.bam; do
    NanoCaller --bam $READS --ref $REF --preset ont --mode snps --mincov 20 --min_allele_freq 0.1 --output $WORKDIR/variant_calling/nanocaller/Polio2 --prefix $(basename $READS).nanocaller
done

for READS in $WORKDIR/data_preprocessing/Polio2/*.sort.bam; do
    NanoCaller --bam $READS --ref $REF2 --preset ont --mode snps --mincov 20 --min_allele_freq 0.1 --output $WORKDIR/variant_calling/nanocaller/Polio2 --prefix $(basename $READS).nanocaller
done

for READS in $WORKDIR/data_preprocessing/Polio3/*.sort.bam; do
    NanoCaller --bam $READS --ref $REF3 --preset ont --mode snps --mincov 20 --min_allele_freq 0.1 --output $WORKDIR/variant_calling/nanocaller/Polio3 --prefix $(basename $READS).nanocaller
done
##################################################

#for READS in $PWD/*.snps.vcf.gz; do
#    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[ %DP]\t%FQ\n' $READS > $(basename $READS).txt
#done

for READS in $PWD/samtools_ampliconclip_both_ends/*.sort.bam; do
    NanoCaller --bam $READS --ref $PWD/Referenz_poliovirus_Sabin_1_converted.fasta --preset ont --mode snps --mincov 20 --min_allele_freq 0.1 --output $PWD/nanocaller_run_clipped_3.0.0 --prefix $(basename $READS)
done
 
for READS in $PWD/samtools_ampliconclip_both_ends/*.sort.bam; do
    NanoCaller --bam $READS --ref $PWD/Referenz_poliovirus_Sabin_1_converted.fasta --preset ont --snp_model ONT-HG002_r10.3 --mode snps --mincov 20 --min_allele_freq 0.1 --output $PWD/nanocaller_run_clipped_3.4.1 --prefix $(basename $READS)
done
```
### format v.3.4.0 vcf such that ALT field is correct
```bash
conda activate formatNanocallerVcf
for READS in $PWD/*.snps.vcf.gz; do
    python3 $WORKDIR/formatNanocallerVcf.py $READS $(basename $READS).corrected.vcf.gz
done

conda activate bcftools
for READS in $PWD/*.snps.vcf.gz; do
    bcftools head $READS > $(basename $READS).header.txt
    bcftools reheader $(basename $READS).corrected.vcf.gz --header $(basename $READS).header.txt > $(basename $READS).final.vcf.gz 
done

for READS in $PWD/*.final.vcf.gz; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%secondary_prob\n' $READS > $(basename $READS ).txt
done

```

## iVar variant calling

```bash
# 2nd run with corrected bed file
wget -O genes.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=AY184219"
for READS in  $WORKDIR/data_preprocessing/Polio1/*.sort.bam; do
    samtools mpileup -aa -A -d 60000 -B -Q 0 $READS | ivar variants -p $(basename $READS) -q 20 -t 0.1 -r $REF -g genes.gff -m 1000
done

wget -O genes.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=AY184220"
for READS in  $WORKDIR/data_preprocessing/Polio2/*.sort.bam; do
    samtools mpileup -aa -A -d 60000 -B -Q 0 $READS | ivar variants -p $(basename $READS) -q 20 -t 0.1 -r $REF2 -g genes.gff -m 1000
done

wget -O genes.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=AY184221"
for READS in  $WORKDIR/data_preprocessing/Polio3/*.sort.bam; do
    samtools mpileup -aa -A -d 60000 -B -Q 0 $READS | ivar variants -p $(basename $READS) -q 20 -t 0.1 -r $REF3 -g genes.gff -m 1000
done

# convert ivar tsv to vcf 
wget "https://github.com/jts/ncov-tools/blob/7a19778c644594fe17ce4fc703560e97b149aa60/workflow/scripts/ivar_variants_to_vcf.py"

for READS in $PWD/*.tsv; do
    python $WORKDIR/variant_calling/ivar/ivar_variants_to_vcf.py $READS $(basename $READS ).ivar.vcf
done

# remove multiallelic SNPs 
for READS in $PWD/*.vcf; do
    bcftools view --max-alleles 2 --exclude-types indels $READS > $(basename $READS ).filtered.vcf
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ALT_FREQ]\n' $(basename $READS ).filtered.vcf > $(basename $READS ).filtered.vcf.txt
done


```

## Clair3 variant calling
```bash
# docker run --rm -it -v $PWD:$PWD -w $PWD hkubal/clair3:latest /bin/bash

MODEL_PATH=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/Nanopore_analysis/clai3_run/models/r1041_e82_400bps_sup_v420
conda activate clair3

for READS in $WORKDIR/data_preprocessing/Polio1/*.sort.bam; do
    mkdir $PWD/$(basename $READS)
    run_clair3.sh \
    --bam_fn=$READS \
    --ref_fn=$REF\
    --threads=2 \
    --platform='ont' \
    --model_path=$MODEL_PATH \
    --chunk_size=2500 \
    --include_all_ctgs \
    --min_coverage=20 \
    --no_phasing_for_fa \
    --haploid_sensitive \
    --snp_min_af=0.1 \
    --call_snp_only \
    --var_pct_full=1 \
    --ref_pct_full=1 \
    --output=$PWD/$(basename $READS) 
done

#### 2nd run with corrected bed file
for READS in $WORKDIR/data_preprocessing/Polio1/*.sort.bam; do
    #mkdir $PWD/$(basename $READS)
    run_clair3.sh \
    --bam_fn=$READS \
    --ref_fn=$REF\
    --threads=2 \
    --platform='ont' \
    --model_path=$MODEL_PATH \
    --chunk_size=2500 \
    --include_all_ctgs \
    --min_coverage=20 \
    --no_phasing_for_fa \
    --haploid_sensitive \
    --snp_min_af=0.1 \
    --call_snp_only \
    --var_pct_full=1 \
    --ref_pct_full=1 \
    --output=$PWD/$(basename $READS) 
done

for READS in $WORKDIR/data_preprocessing/Polio2/*.sort.bam; do
    mkdir $PWD/$(basename $READS)
    run_clair3.sh \
    --bam_fn=$READS \
    --ref_fn=$REF2\
    --threads=2 \
    --platform='ont' \
    --model_path=$MODEL_PATH \
    --chunk_size=2500 \
    --include_all_ctgs \
    --min_coverage=20 \
    --no_phasing_for_fa \
    --haploid_sensitive \
    --snp_min_af=0.1 \
    --call_snp_only \
    --var_pct_full=1 \
    --ref_pct_full=1 \
    --output=$PWD/$(basename $READS) 
done

for READS in $WORKDIR/data_preprocessing/Polio3/*.sort.bam; do
    mkdir $PWD/$(basename $READS)
    run_clair3.sh \
    --bam_fn=$READS \
    --ref_fn=$REF3\
    --threads=2 \
    --platform='ont' \
    --model_path=$MODEL_PATH \
    --chunk_size=2500 \
    --include_all_ctgs \
    --min_coverage=20 \
    --no_phasing_for_fa \
    --haploid_sensitive \
    --snp_min_af=0.1 \
    --call_snp_only \
    --var_pct_full=1 \
    --ref_pct_full=1 \
    --output=$PWD/$(basename $READS) 
done


find ./  -name "merge_output.vcf.gz" > vcf_files.txt
cat vcf_files.txt | while read line; do 
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[ %DP]\t[ %AF]\n' $line >  $(dirname $line)/$(basename $line).txt 
done

```

## create .tabular files for heatmap
```bash
for READS in $PWD/*.annotation.covered.af.vcf.gz; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%DP\t%AF\n' $READS > $(basename $READS).tabular
    tr '|' '\t' < $(basename $READS).tabular > $(basename $READS).1.tabular
done
# covpipe *.annotation.covered.af.vcf.gz
# nanopore *.ann.vcf


# Prepare clair3 for HEATMAP 
find -name "merge_output.vcf.gz" > vcf_files.txt
while read line ; do
    vcf_name=$(echo $(basename $(dirname $line)) | cut -d'.' -f1)
    tab=".tabular"
    tab2="_clair3"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[ %DP]\t[ %AF]\n' $line > "$vcf_name$tab"
    tr '|' '\t' < "$vcf_name$tab" > "$vcf_name$tab2$tab"
    sed -i '1s/^/CHROM\tPOS\tREF\tALT\tDP\tAF\n/' "$vcf_name$tab2$tab"
done < $PWD/vcf_files.txt

# Prepare iVAR for HEATMAP 
for line in $PWD/*.filtered.vcf; do
    vcf_name=$(echo $(basename $line) | cut -d'.' -f1)
    tab=".tabular"
    tab2="_iVAR"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t[%ALT_FREQ]\n' $line > "$vcf_name$tab"
    tr '|' '\t' < "$vcf_name$tab" > "$vcf_name$tab2$tab"
    sed -i '1s/^/CHROM\tPOS\tREF\tALT\tDP\tAF\n/' "$vcf_name$tab2$tab"
done 

# Prepare medaka for HEATMAP 
for line in $PWD/*snp.vcf; do
    vcf_name=$(echo $(basename $line) | cut -d'.' -f1)
    tab=".tabular"
    tab2="_medaka_snp"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[ %GQ]\t%primary_prob\n' $line > "$vcf_name$tab"
    tr '|' '\t' < "$vcf_name$tab" > "$vcf_name$tab2$tab"
    sed -i '1s/^/CHROM\tPOS\tREF\tALT\tDP\tAF\n/' "$vcf_name$tab2$tab"
done 

# Prepare nanocaller for HEATMAP 
for line in $PWD/*snps.vcf.gz; do
    vcf_name=$(echo $(basename $line) | cut -d'.' -f1)
    tab=".tabular"
    tab2="_nanocaller"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[ %DP]\t[ %FQ]\n' $line > "$vcf_name$tab"
    tr '|' '\t' < "$vcf_name$tab" > "$vcf_name$tab2$tab"
    sed -i '1s/^/CHROM\tPOS\tREF\tALT\tDP\tAF\n/' "$vcf_name$tab2$tab"
done 

# Prepare piranha-medaka_haploid_variant for HEATMAP -----> Whole Genome run with our own REF123.fasta
#### WARNING -> for plotting purposes GQ and GT are taken as place holder for DP and AF.
#### Medaka does not output DP and AF !!!
for line in $PWD/barcode*/variant_calls/*.vcf; do
    vcf_name=$(echo  $(dirname $line) | cut -d'/' -f9)
    tab=".tabular"
    tab2="_piranha_medaka_haploid_variant_"
    tab3=$(basename $line)

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[ %GQ]\t[ %GT]\n' $line > "$vcf_name$tab"
    tr '|' '\t' < "$vcf_name$tab" > "$vcf_name$tab2$tab3$tab"
    sed -i '1s/^/CHROM\tPOS\tREF\tALT\tDP\tAF\n/' "$vcf_name$tab2$tab3$tab" 
done 

# Prepare Covpipe2-freebayes for HEATMAP
for line in $PWD/results/03-Variant-Calling/SiB*/*.filtered.gt_adjust.filtered_indels.vcf.gz; do
    vcf_name=$(echo  $(basename $line) | cut -d'.' -f1)
    tab=".tabular"
    tab2="_Illumina_freebayes"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AF\n' $line > "$vcf_name$tab"
    tr '|' '\t' < "$vcf_name$tab" > "$vcf_name$tab2$tab"
    sed -i '1s/^/CHROM\tPOS\tREF\tALT\tDP\tAF\n/' "$vcf_name$tab2$tab"
done 

# Prepare Sanger vcf for HEATMAP 
#### WARNING -> for plotting purposes DP and AF have placeholder values.
for line in $PWD/*.fasta.vcf; do
    vcf_name=$(echo $(basename $line) | cut -d'.' -f1)
    tab=".tabular"
    tab2="sanger"
    tab3="_"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $line > "$vcf_name$tab"
    tr '|' '\t' < "$vcf_name$tab" > "$vcf_name$tab3$tab"
    awk -F'\t' -vOFS='\t' '{ $6 = "1";  $5 = "100"; $1 = "AY184221" }1' < "$vcf_name$tab3$tab" > "$vcf_name$tab3$tab2$tab"
    sed -i '1s/^/CHROM\tPOS\tREF\tALT\tDP\tAF\n/' "$vcf_name$tab3$tab2$tab"
    rm "$vcf_name$tab" "$vcf_name$tab3$tab"
done 
```

## Run Piranha

```bash
FASTQPIRANHA=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/data/2023-08-18-PolioPrimer-ONT-runID114/piranha_demultiplexed

piranha -i $FASTQPIRANHA -b $PWD/barcodes.csv -r $REF -m wg --medaka-model r1041_e82_400bps_sup_v4.0.0 --min-read-length 600 --max-read-length 1000 -s environmental --min-read-depth 30 --save-config

piranha -i $FASTQPIRANHA -b $PWD/barcodes.csv --analysis-mode wg --medaka-model r1041_e82_400bps_sup_v4.0.0 --min-read-length 600 --max-read-length 1000 -s environmental --no-temp --output-prefix analysis_with_phylo --save-config

### Piranha run 08-11-2023
piranha -i $FASTQPIRANHA -b $PWD/barcodes.csv -r $REF123 -m wg --medaka-model r1041_e82_400bps_sup_v4.0.0 --min-read-length 600 --max-read-length 1000 -s environmental --min-read-depth 30 --no-temp --output-prefix analysis_with_ref123 --save-config
```

## rename barcode files from Piranha
```bash
while IFS=, read col1 col2 ; do
    tab_name=$(basename $col1"_piranha_medaka_haploid_variant_"*)
    echo $tab_name
    sib=$col2
    echo "${sib}_piranha_medaka_haploid_variant.tabular"
    #mv $tab_name $col2""
    #rename 's/^barcode01/SiB155/' barcode19*
done < barcodes.csv
```

# De Novo Assembly
## Run Flye on read length filtered data
```bash
conda activate flye

for READS in $WORKDIR/fastq_barcode_merged/run_polio3/*.filtered.fastq.gz; do
    mkdir $PWD/$(basename $READS)
    flye --nano-raw $READS -o $PWD/$(basename $READS) -t 8  --meta --min-overlap 1000 -g 7400
done

for READS in $WORKDIR/data_preprocessing/Polio3/*.filtered.fastq.gz; do
    mkdir $PWD/$(basename $READS)
    flye --nano-raw $READS -o $PWD/$(basename $READS) -t 8  --meta --min-overlap 1000 -g 7400
done
## get contig and contig sizes
conda activate bioawk
for READS in $PWD/20230403_1517_X5_FAU91805_e08a3e32_AS-28*/*contigs.fasta; do
    echo $READS 
    bioawk -c fastx '{ print $name, length($seq) }' < $READS 
done > contig_sizes.txt

## try de novo assembly with spades (made for illumina and not nanopore)
spades.py --rnaviral -o polio2 -s $WORKDIR/fastq_barcode_merged/run_polio2/*.filtered.fastq.gz
spades.py --rnaviral -o polio1 -s $WORKDIR/fastq_barcode_merged/run_polio1/20230403_1517_X5_FAU91805_e08a3e32_AS-2858_fastq_pass.filtered.fastq.gz

# de novo assembly with canu
conda activate canu

for READS in $WORKDIR/fastq_barcode_merged/run_polio3/*.filtered.fastq.gz; do
    canu -assemble -nanopore $READS -d $PWD/$(basename $READS) -p $(basename $READS) genomeSize=7.4k -corrected
done
```
# Check for soft-clipped alignments
```bash
SAMDIR=$WORKDIR/fastq_barcode_merged/run_polio1
for SAM in $SAMDIR/*.sam; do
    python samclipy.py < $SAM > $SAMDIR/$(basename $SAM).softclip.filtered.sam 
done > $SAMDIR/samclipy.log

for READS in $SAMDIR/*.softclip.filtered.sam; do
    samtools view -bS $READS | samtools sort - -@ 2 -o $SAMDIR/$(basename $READS).sorted.bam
    samtools index $SAMDIR/$(basename $READS).sorted.bam
done

# get reads with H (hard-) or S (soft-) clip in region 5400-6000 ---> coverage drop in PV1
samtools view -h 20230403_1517_X5_FAU91805_e08a3e32_AS-2864_fastq_pass.filtered.fastq.gz.sorted.bam "AY184219:5400-6000" | awk '$6 ~ /H/ || $1 ~ /^@/' | samtools view -bS - > 20230403_1517_X5_FAU91805_e08a3e32_AS-2864_fastq_pass.filtered.fastq.gz.sorted.bam_CIGAR_H-output.bam

# get number of reads in bam/sam file
samtools view -c 20230403_1517_X5_FAU91805_e08a3e32_AS-2871_fastq_pass.filtered.fastq.gz.sam
# samtools view -h input.bam | awk '$6 !~ /S/ || $1 ~ /^@/' | samtools view -bS - > not-n-output.bam
```

# create consensus sequence -> did not work
```bash
# make primer fa file
PRIMER_DIR=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/primer/polio2
while read line ; do
    set $line 
    echo ">$2"
    echo ${4^^} 
done < <(tail -n +2 $PWD/0_aligned_consensus_masked.fasta_900.primer.tsv) > primer.fa

# run consensus
gunzip $WORKDIR/data_preprocessing/Polio2/*.filtered.fastq.gz 
for READS in $WORKDIR/data_preprocessing/Polio2/*.filtered.fastq; do
    echo $READS 
    NGSpeciesID --ont --consensus --racon --racon_iter 3 --fastq $READS --outfolder $PWD --primer_file $PRIMER_DIR/primer.fa
done
 
for file in *.fastq; do
bn=`basename $file .fastq`
NGSpeciesID --ont --consensus --sample_size 500 --m 800 --s 100 --medaka --primer_file primers.txt --fastq $file --outfolder ${bn}
done
gzip $WORKDIR/data_preprocessing/Polio2/*.filtered.fastq 


```