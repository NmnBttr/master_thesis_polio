# Short analysis for Polio SANGER Data

## SETUP
```bash
srun --cpus-per-task=4 --mem=50GB --gres=local:30 --gpus=1 --pty --time=08:00:00 bash -i 

WORKDIR=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/Sanger_analysis
FASTA=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/SANGER_VP1_sequences
REF=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/References/Referenz_Poliovirus_Sabin_1.fasta
REF2=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/References/Referenz_Poliovirus_Sabin_2.fasta
REF3=/scratch/projekte/MF1_Bioinformatics/poliovirus-primer-evaluation/References/Referenz_Poliovirus_Sabin_3.fasta
```
# prepare fasta file
```bash
# split multi fasta into individual samples
cat Polio_Sequenzen_Tiled_PCR.fasta | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($1,2,6) "_VP1.fasta")}
        print $0 >> filename
        close(filename)
}'

cat 3D_Sanger.fasta  | awk '{
        if (substr($0, 1, 1)==">") {filename=(substr($1,2,6) "_3D.fasta")}
        print $0 >> filename
        close(filename)
}'
```

# MSA (pairwise alignment to reference genomes) using MAFFT
```bash
conda activate mafft

# mafft on individual sample fasta files -> pairwise alignment to Reference genome
cd Polio1
for fasta in $FASTA/Polio1/*.fasta; do
  mafft --auto --quiet --addfragments $fasta --reorder --kimura 1 --thread -1 $REF > $(basename $fasta).aligned.fasta
done
cd ..
cd Polio2
for fasta in $FASTA/Polio2/*.fasta; do
  mafft --auto --quiet --addfragments $fasta --reorder --kimura 1 --thread -1 $REF2 > $(basename $fasta).aligned.fasta
done
cd ..
cd Polio3
for fasta in $FASTA/Polio3/*.fasta; do
  mafft --auto --quiet --addfragments $fasta --reorder --kimura 1 --thread -1 $REF3 > $(basename $fasta).aligned.fasta
done
cd ..
```

# MSA Fasta to VCF https://github.com/sanger-pathogens/snp-sites/tree/52c98cb3e0ed0d336b24b27a5c0f3da4cbe44e71

```bash
conda activate snp-sites
# per sample
for FASTA in $WORKDIR/MSA_mafft/Polio3/*.fasta; do
    snp-sites -cv -o $(basename $FASTA).vcf $FASTA
done
conda deactivate

conda activate sambcfenv
# splitting the main vcf into individual sample vcf   
# not for individual samples 
for file in *.vcf; do
  for sample in `bcftools query -l $file`; do
    bcftools view -c1 -Oz -s $sample -o ${file/.vcf/.$sample.vcf} $file
  done
done 
```
