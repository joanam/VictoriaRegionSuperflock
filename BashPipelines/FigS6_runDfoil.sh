#!/bin/bash

prefix=$1
i=$2
P1=$3
P2=$4
P3=$5
P4=$6
P5=$7

module load python/2.7

file=allGenomes.1kbThinned.lakes.tree.groups
file=allGenomes.LDpruned.lakes.tree.groups
file=allGenomes.1kbThinned.lakes.tree.groups.withSRA.chr1-21

vcftools --gzvcf /cluster/work/gdc/shared/p500/victoriaGenomes/LVRS/newLVRS/vcf/$file.vcf.gz --indv $P1 --indv $P2 --indv $P3 --indv $P4 --indv $P5 \
  --recode --max-missing 1 --mac 1 --stdout > $prefix.$i.vcf

# Convert to phylip alignment with my python script
python2.7 `which vcf2phylip.py` -i $prefix.$i.vcf -o $prefix.$i.phy

# Convert phylip to fasta
awk '{if(NR>1) print ">"$1"\n"$2}' $prefix.$i.phy  > $prefix.$i.fa

# Compute Dfoil pattern counts
fasta2dfoil.py  --out $prefix.$i.dfoil.counts --names $P1,$P2,$P3,$P4,$P5 $prefix.$i.fa

# run Dfoil
python `which dfoil.py` --infile $prefix.$i.dfoil.counts --out $prefix.$i.dfoil --mode dfoilalt

# Delete intermediary files
rm $prefix.$i.phy $prefix.$i.fa $prefix.$i.vcf
