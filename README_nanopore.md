# SETUP

```bash
WORKDIR=/home/namuun/Documents/Polio_Nanopore_Analysis
FASTQ=/home/namuun/Documents/CoVpipe_Polio/data/230215_GI1-3_Run23-044/kreibichj
REF=/home/namuun/Documents/Polio_Nanopore_Analysis/Referenz_poliovirus_Sabin_1_converted.fasta
PRIMER=/home/namuun/Downloads/Poliovirus_1/results/polio/final_primers/0_aligned_consensus_masked/primer_polioS1.bed
GFF=$WORKDIR/genes.gff
# 0_aligned_consensus_masked.fasta_1000.primer.bed
# primer_polioS1.bed

KRAKENDB=/home/namuun/Documents/Polio_Nanopore_Analysis/minikraken_20141208/database.kdb.tar.gz
```

# prepare reads, 2 samples each 2 seqs
```bash
for READS in $FASTQ/barcode*; do BN=$(basename $READS); cat $READS/*.fastq.gz > $BN.fastq.gz; done

cat barcode21_SiB_151_23-01715.fastq.gz barcode22_SiB_152_23-01716.fastq.gz > SiB_151_152_23-0716.fastq.gz
cat barcode23_SiB_153_23-01717.fastq.gz barcode24_SiB_154_23-01718.fastq.gz > SiB_153_154_23-01717.fastq.gz
```

# Filter reads for the expected amplicon length (+/- 100 bp)

```bash
cd Filtlong

for READS in $WORKDIR/*.fastq.gz; do
    bin/filtlong -a $REF --min_length 100 --keep_percent 90 --trim $READS | gzip > $(basename $READS .fastq.gz).filtered.fastq.gz
done

# bin/filtlong -a $REF --min_length 100 --keep_percent 90 --trim /home/namuun/Documents/Polio_Nanopore_Analysis/SiB_151_152_23-0716.fastq.gz | gzip > SiB_151_152_23-0716_filtered.fastq.gz

# bin/filtlong -a $REF --min_length 100 --keep_percent 90 --trim /home/namuun/Documents/Polio_Nanopore_Analysis/SiB_151_152_23-0716.fastq.gz | gzip > SiB_151_152_23-0716_filtered.fastq.gz
```

# Map the filtered reads to a reference genome
```bash
conda activate minimap2

for READS in /home/namuun/Documents/Polio_Nanopore_Analysis/*.fastq.gz; do
    minimap2 -ax map-ont $REF $READS > $(basename $READS .fastq.gz).sam
    samtools view -bS $(basename $READS .fastq.gz).sam | samtools sort - -@ 2 -o $(basename $READS .fastq.gz).sorted.bam
    samtools index $(basename $READS .fastq.gz).sorted.bam
done

```

# Primer clippping

```bash
for BAMF in $WORKDIR/*.sorted.bam; do
    samtools ampliconclip -b $PRIMER --strand --both-ends $BAMF -o $(basename $BAMF).clipped.bam | tee $(basename $BAMF).samtools.log
    samtools sort $(basename $BAMF).clipped.bam > $(basename $BAMF ).clipped.sort.bam
    samtools index $(basename $BAMF ).clipped.sort.bam
    samtools fastq $(basename $BAMF ).clipped.sort.bam > $(basename $BAMF ).clipped.sort.bam.fastq
done
    
# bamclipper.sh -b /home/namuun/Documents/Polio_Nanopore_Analysis/SiB_151_152_23-0716_filtered.sorted.bam -p $PRIMER
```

# check primer clipping
```bash
for BAMF in $PWD/*.fastq; do
    while read line ; do
        set $line 
        grep ${4^^}$ $BAMF
        grep ^${4^^} $BAMF
    done < <(tail -n +2 /home/namuun/Downloads/Poliovirus_1/results/polio/final_primers/0_aligned_consensus_masked/0_aligned_consensus_masked.fasta_1000.primer.tsv) > $(basename $BAMF).txt
done
```


# Variant calling with basecall_model_version_id=dna_r1041_e82_260bps_sup@v3.5.2
# Run Medaka
# ATTENTION: it is always good to assign an appropriate Medaka model based on the performed basecalling! 
```bash

for READS in $WORKDIR/samtools_ampliconclip_both_ends/*.sort.bam; do
    medaka consensus $READS --model r1041_e82_260bps_sup_v4.0.0 --threads 2 $(basename $READS .output_consensus).hdf;
    medaka snp $REF $(basename $READS .output_consensus).hdf $(basename $READS ).snp.vcf --verbose --threshold 0.1;       
done

for READS in $WORKDIR/samtools_ampliconclip_both_ends/*.sort.bam; do
    medaka consensus $READS --model r1041_e82_260bps_sup_g632 --threads 2 $(basename $READS .output_consensus).hdf;
    medaka snp $REF $(basename $READS .output_consensus).hdf $(basename $READS ).snp.vcf --verbose --threshold 0.1;       
done

for READS in $WORKDIR/samtools_ampliconclip_both_ends/*.sort.bam; do
    medaka consensus $READS --model r1041_e82_260bps_sup_v4.1.0 --threads 2 $(basename $READS .output_consensus).hdf;
    medaka snp $REF $(basename $READS .output_consensus).hdf $(basename $READS ).snp.vcf --verbose --threshold 0.1;       
done

for READS in $WORKDIR/samtools_ampliconclip_both_ends/*.sort.bam; do
    medaka consensus $READS --model r1041_e82_260bps_sup_variant_v4.1.0 --threads 2 $(basename $READS .output_consensus).hdf;
    medaka snp $REF $(basename $READS .output_consensus).hdf $(basename $READS ).snp.vcf --verbose --threshold 0.1;       
done

for READS in $WORKDIR/samtools_ampliconclip_both_ends/*.sort.bam; do
    medaka consensus $READS --model r1041_e82_260bps_sup_variant_g632 --threads 2 $(basename $READS .output_consensus).hdf;
    medaka snp $REF $(basename $READS .output_consensus).hdf $(basename $READS ).snp.vcf --verbose --threshold 0.1;       
done

# cat vcf_files.txt | while read line; do bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%secondary_prob\n' $line  > vcf_files.txt; done 

for READS in $PWD/*.snp.vcf; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%secondary_prob\n' $READS > $(basename $READS ).txt
done

#medaka variant $REF $(basename $READS .output_consensus).hdf $(basename $READS ).variant.vcf; 
# -m r1041_e82_260bps_sup_g632 --> doesnt work well
# -m r1041_e82_260bps_sup_v4.0.0 --> works best
# r1041_e82_260bps_sup_v4.1.0 --> works okay
# r1041_e82_260bps_sup_variant_g632 --> doesnt work well
# r1041_e82_260bps_sup_variant_v4.1.0 --> doesnt work well
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
sudo chmod 666 /var/run/docker.sock

docker run --rm -it -v $PWD:$PWD -w $PWD genomicslab/nanocaller:3.4.1 /bin/bash

for READS in $PWD/samtools_ampliconclip_both_ends/*.sort.bam; do
    NanoCaller --bam $READS --ref $PWD/Referenz_poliovirus_Sabin_1_converted.fasta --preset ont --mode snps --mincov 20 --min_allele_freq 0.1 --output $PWD/nanocaller_run_clipped_3.4.0 --prefix $(basename $READS)
done

for READS in $PWD/*.snps.vcf.gz; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FQ\n' $READS > $(basename $READS ).txt
done

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
for READS in $WORKDIR/samtools_ampliconclip_both_ends/*.sort.bam; do
    samtools mpileup -aa -A -d 60000 -B -Q 0 $READS | ivar variants -p $(basename $READS) -q 20 -t 0.1 -r $REF -g $GFF -m 1000
done

# convert ivar tsv to vcf 
wget "https://github.com/jts/ncov-tools/blob/7a19778c644594fe17ce4fc703560e97b149aa60/workflow/scripts/ivar_variants_to_vcf.py"
python ivar_variants_to_vcf.py /home/namuun/Documents/Polio_Nanopore_Analysis/ivar_run_clipped/SiB_151_152_23-0716_filtered.sorted.bam.clipped.sort.bam.tsv /home/namuun/Documents/Polio_Nanopore_Analysis/ivar_run_clipped/SiB_151_152_23-0716_filtered.sorted.bam.clipped.sort.bam.vcf


for READS in $PWD/*.tsv; do
    python ivar_variants_to_vcf.py $READS > $(basename $READS ).vcf
done
for READS in $PWD/*filtered.vcf; do
    
done

# remove multiallelic SNPs 
for READS in $PWD/*.sort.bam.vcf; do
    bcftools view --max-alleles 2 --exclude-types indels $READS > $(basename $READS ).filtered.vcf
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%ALT_FREQ]\n' $(basename $READS ).filtered.vcf > $(basename $READS ).filtered.vcf.txt
done


```

## Clair3 variant calling
```bash
docker run --rm -it -v $PWD:$PWD -w $PWD hkubal/clair3:latest /bin/bash

for READS in $PWD/samtools_ampliconclip_both_ends/*.sort.bam; do
    mkdir $PWD/clair3_run_clipped/$(basename $READS)_run16

    /opt/bin/run_clair3.sh \
    --bam_fn=$READS \
    --ref_fn=/home/namuun/Documents/Polio_Nanopore_Analysis/Referenz_poliovirus_Sabin_1_converted.fasta\
    --threads=2 \
    --platform='ont' \
    --model_path="/home/namuun/Documents/Polio_Nanopore_Analysis/rerio_clair3_models/rerio/clair3_models/r1041_e82_260bps_sup_v410" \
    --chunk_size=2500 \
    --include_all_ctgs \
    --min_coverage=20 \
    --no_phasing_for_fa \
    --haploid_sensitive \
    --snp_min_af=0.1 \
    --call_snp_only \
    --output=/home/namuun/Documents/Polio_Nanopore_Analysis/clair3_run_clipped/$(basename $READS)_run16 
done
```
--min_coverage=20 --call_snp_only --indel_min_af=0.1 --haploid_precise
--chunk_size=1000 --include_all_ctgs --no_phasing_for_fa
# model r941_prom_sup_g5014: empty vcf
# model r941_prom_hac_g360+g422: empty vcf
# model r941_prom_hac_g238: empty vcf

/opt/bin/scripts/clair3_c_impl.sh --bam_fn /home/namuun/Documents/Polio_Nanopore_Analysis/samtools_ampliconclip_both_ends/SiB_151_152_23-0716_filtered.sorted.bam.clipped.sort.bam --ref_fn /home/namuun/Documents/Polio_Nanopore_Analysis/Referenz_poliovirus_Sabin_1_converted.fasta --threads 2 --model_path /home/namuun/Documents/Polio_Nanopore_Analysis/rerio_clair3_models/rerio/clair3_models/r1041_e82_260bps_sup_v410 --platform ont --output /home/namuun/Documents/Polio_Nanopore_Analysis/clair3_run_clipped/SiB_151_152_23-0716_filtered.sorted.bam.clipped.sort.bam_run5 --bed_fn=EMPTY --vcf_fn=EMPTY --ctg_name=EMPTY --sample_name=SAMPLE --chunk_num=0 --chunk_size=1000 --samtools=samtools --python=python3 --pypy=pypy3 --parallel=parallel --whatshap=whatshap --qual=2 --var_pct_full=0.7 --ref_pct_full=0.1 --var_pct_phasing=0.7 --snp_min_af=0.08 --indel_min_af=0.15 --min_mq=5 --min_coverage=20 --min_contig_size=0 --pileup_only=False --gvcf=False --fast_mode=False --call_snp_only=True --print_ref_calls=False --haploid_precise=False --haploid_sensitive=False --include_all_ctgs=True --no_phasing_for_fa=False --pileup_model_prefix=pileup --fa_model_prefix=full_alignment --remove_intermediate_dir=False --enable_phasing=False --enable_long_indel=False --keep_iupac_bases=False --use_gpu=False --longphase_for_phasing=False --longphase=EMPTY --use_whatshap_for_intermediate_phasing=True --use_longphase_for_intermediate_phasing=False --use_whatshap_for_final_output_phasing=False --use_longphase_for_final_output_phasing=False --use_whatshap_for_final_output_haplotagging=False

### flye --nano-raw SiB_151_152_23-0716_filtered.fastq -o flye_output -t 8 --meta -g 7k

# Annotation with SnpEff

```bash
# retrieve genebank file for Polio 1
# wget -O genes.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=AY184219"
# gzip genes.gff
# esearch -db nucleotide -query "AY184219.1" | efetch -format fasta > AY184219.1.fa

# Building a database from GenBank files https://pcingola.github.io/SnpEff/se_build_db/#step-1-configure-a-new-genome
#  https://www.ncbi.nlm.nih.gov/nuccore/AY184219.1
# build database for SnpEff
java -jar snpEff.jar build -genbank -v AY184219.1

for READS in $PWD/*.snps.vcf; do
    sed -i 's/AY184219/AY184219.1/' $READS
done

#freebayes
mkdir SiB_151_152
mkdir SiB_153_154
cd SiB_15....
java -Xmx8g -jar /home/namuun/Documents/Polio_Nanopore_Analysis/snpEff/snpEff.jar AY184219.1 -formatEff $WORKDIR/freebayes_run_clipped/SiB_151_152_23-0716_filtered.sorted.bam.clipped.sort.bam.vcf > SiB_151_152.ann.vcf

java -Xmx8g -jar /home/namuun/snpEff/snpEff.jar AY184219.1 -formatEff $WORKDIR/freebayes_run_clipped/SiB_153_154_23-01717_filtered.sorted.bam.clipped.sort.bam.vcf > SiB_153_154.ann.vcf

# medaka
java -Xmx8g -jar /home/namuun/snpEff/snpEff.jar AY184219.1 -formatEff $WORKDIR/freebayes_run_clipped/SiB_153_154_23-01717_filtered.sorted.bam.clipped.sort.bam.snp.vcf > SiB_153_154.ann.vcf

# nanocaller 
gunzip *final.vcf.gz
for READS in $PWD/*.final.vcf; do
    sed -i 's/AY184219/AY184219.1/' $READS
done
## 3.4.0
java -Xmx8g -jar /home/namuun/snpEff/snpEff.jar AY184219.1 -formatEff $WORKDIR/nanocaller_run_clipped_3.4.0_backup/SiB_153_154_23-01717_filtered.sorted.bam.clipped.sort.bam.unfiltered.snps.vcf.gz.final.vcf > SiB_153_154.ann.vcf
## 3.0.0 
java -Xmx8g -jar /home/namuun/snpEff/snpEff.jar AY184219.1 -formatEff $WORKDIR/nanocaller_run_clipped_3.0.0_copy/SiB_153_154_23-01717_filtered.sorted.bam.clipped.sort.bam.snps.vcf > SiB_153_154.ann.vcf

# Illumina
java -Xmx8g -jar /home/namuun/snpEff/snpEff.jar AY184219.1 -formatEff /home/namuun/Documents/CoVpipe_Polio/Illumina_complete_run3/results/03-Variant-Calling/SiB_147_148_S_6_7/SiB_147_148_S_6_7.vcf > SiB_151_152_Illumina_freebayes.ann.vcf  

java -Xmx8g -jar /home/namuun/snpEff/snpEff.jar AY184219.1 -formatEff /home/namuun/Documents/CoVpipe_Polio/Illumina_complete_run3/results/03-Variant-Calling/SiB_149_150_S_8_9/SiB_149_150_S_8_9.vcf > SiB_153_154_Illumina_freebayes.ann.vcf

# ivar
java -Xmx8g -jar /home/namuun/Documents/Polio_Nanopore_Analysis/snpEff/snpEff.jar AY184219.1 -formatEff /home/namuun/Documents/Polio_Nanopore_Analysis/ivar_run_clipped/SiB_151_152_23-0716_filtered.sorted.bam.clipped.sort.bam.vcf.filtered.vcf > SiB_151_152_ivar.ann.vcf  

java -Xmx8g -jar /home/namuun/Documents/Polio_Nanopore_Analysis/snpEff/snpEff.jar AY184219.1 -formatEff /home/namuun/Documents/Polio_Nanopore_Analysis/ivar_run_clipped/SiB_153_154_23-01717_filtered.sorted.bam.clipped.sort.bam.vcf.filtered.vcf > SiB_153_154_ivar.ann.vcf

```
# medaka
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%secondary_prob\t%EFF\n' 
# freebayes
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%DP\t%AB\t%EFF\n'
# nanocaller 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[ %DP]\t%FQ\t%EFF\n'
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[ %DP]\t[ %FQ]\t%EFF\n' SiB_151_152.ann.vcf > SiB_151_152.ann.tabular
# ivar
for READS in $PWD/*.ann.vcf; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%DP\t[%ALT_FREQ]\t%EFF\n' $READS > $(basename $READS).tabular
    tr '|' '\t' < $(basename $READS).tabular > $(basename $READS).1.tabular
done

for READS in $PWD/*.ann.vcf; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t[%ALT_FREQ]\n' $READS
done


##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype [ | ERRORS | WARNINGS ] )' ">

CHROM	POS	REF	ALT	FILTER	DP	AF	EFF[*].EFFECT	EFF[*].IMPACT	EFF[*].FUNCLASS	EFF[*].AA	EFF[*].GENE

find -name "*.tabular" > snpeff.tabular.txt

while read line ; do
    tr '|' '\t' < $WORKDIR/$line > $(basename $line).tab.tabular
done < $WORKDIR/snpeff.tabular.txt

for READS in $PWD/*.tabular; do
    sed -i '1s/^/CHROM  POS REF ALT FILTER  DP  AF  EFF[*].EFFECT   EFF[*].IMPACT   EFF[*].FUNCLASS EFF[*].AA   EFF[*].GENE\n/' $READS
done
<!-- java -jar snpEff.jar build -genbank -v AY184219.1
00:00:00 SnpEff version SnpEff 5.1f (build 2023-07-21 13:58), by Pablo Cingolani
00:00:00 Command: 'build'
00:00:00 Building database for 'AY184219.1'
00:00:00 Reading configuration file 'snpEff.config'. Genome: 'AY184219.1'
00:00:00 Reading config file: /home/namuun/snpEff/snpEff.config
00:00:00 done
00:00:00 Chromosome: 'AY184219.1'       length: 7441
00:00:00 Create exons from CDS (if needed): 
............00:00:00 Exons created for 12 transcripts.
00:00:00 Deleting redundant exons (if needed): 
00:00:00        Total transcripts with deleted exons: 0
00:00:00 Collapsing zero length introns (if needed): 
00:00:00        Total collapsed transcripts: 0
00:00:00                Adding genomic sequences to genes: 
00:00:00        Done (1 sequences added).
00:00:00                Adding genomic sequences to exons: 
00:00:00        Done (12 sequences added, 0 ignored).
00:00:00 Finishing up genome
00:00:00 Adjusting transcripts: 
00:00:00 Adjusting genes: 
00:00:00 Adjusting chromosomes lengths: 
00:00:00 Ranking exons: 
00:00:00 Create UTRs from CDS (if needed): 
00:00:00 Remove empty chromosomes: 
00:00:00 Marking as 'coding' from CDS information: 
00:00:00 Done: 0 transcripts marked
00:00:00 
00:00:00 Caracterizing exons by splicing (stage 1) : 

00:00:00 Caracterizing exons by splicing (stage 2) : 
        00:00:00 done.
00:00:00 [Optional] Rare amino acid annotations
WARNING_FILE_NOT_FOUND: Rare Amino Acid analysis: Cannot read protein sequence file '/home/namuun/snpEff/./data/AY184219.1/protein.fa', nothing done.
00:00:00 CDS check: GenBank file format, skipping

00:00:00 Protein check file: '/home/namuun/snpEff/./data/AY184219.1/genes.gbk'

00:00:00 Checking database using protein sequences
00:00:00 Comparing Proteins...
        Labels:
                '+' : OK
                '.' : Missing
                '*' : Error
        ......+.....00:00:00 

        Protein check:  AY184219.1      OK: 1   Not found: 11   Errors: 0       Error percentage: 0.0%
00:00:00 Saving database
00:00:00 Saving sequences for small chromosmes to file '/home/namuun/snpEff/./data/AY184219.1/sequence.bin'
00:00:00 [Optional] Reading regulation elements: GFF
WARNING_FILE_NOT_FOUND: Cannot read optional regulation file '/home/namuun/snpEff/./data/AY184219.1/regulation.gff', nothing done.
00:00:00 [Optional] Reading regulation elements: BED 
00:00:00 Cannot find optional regulation dir '/home/namuun/snpEff/./data/AY184219.1/regulation.bed/', nothing done.
00:00:00 [Optional] Reading motifs: GFF
WARNING_FILE_NOT_FOUND: Cannot open PWMs file /home/namuun/snpEff/./data/AY184219.1/pwms.bin. Nothing done
00:00:00 Done
00:00:00 Logging
00:00:01 Checking for updates...
00:00:02 Done.
 -->

 https://github.com/pcingola/SnpEff/issues/479
https://github.com/pcingola/SnpEff/issues/455
 https://github.com/pcingola/SnpEff/issues/474


## piranha

```bash
piranha -i $FASTQ -b $PWD/barcodes.csv -r $REF -m wg --medaka-model r1041_e82_260bps_sup_v4.0.0
piranha -i $FASTQ -b $PWD/barcodes.csv -m wg --medaka-model r1041_e82_260bps_sup_v4.0.0 --save-config

```