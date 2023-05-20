#!/bin/bash
# Run fineSTRUCTURE


# Phase the dataset (with runBeagle.lsf)
path=/cluster/work/gdc/shared/p500/victoriaGenomes/LVRS/newLVRS/vcf/
file=allGenomes.withSRA.chr1-22.max0.5N.SNPs.forFbranch.100bp

# Generate a plink map file (only once)
addRecombRates.r $path$file.vcf.gz /cluster/work/gdc/shared/p500/ref-genome/puncross.gapsEstimated.lifted.pruned.bed

CHR=chr${LSB_JOBINDEX}

# Generate the Beagle and IBD files
java -Xmx100G -Djava.io.tmpdir=./tempBeagle.2 \
   -jar /cluster/apps/gdc/beagle/4.1/beagle.03May16.862.jar \
   overlap=15000 window=37500 map=$path$file.plink.map \
   gt=$path$file.vcf.gz nthreads=4 niterations=10 ibd=true \
   out=$file.phased.$CHR ibdtrim=1125 chrom=$CHR

# Fix the header (add in chr info)
bsub "cat header.contigs.info <(zgrep -v "^##" allGenomes.withSRA.chr1-22.max0.5N.SNPs.forFbranch.100bp.phased.chr1-6.vcf.gz) | gzip > allGenomes.withSRA.chr1-22.max0.5N.SNPs.forFbranch.100bp.phased.chr1-6.headerFixed.vcf.gz"

# fix headers
for CHR in chr{1..22}
do
 cat header.contigs.info <(zgrep -v '^##' allGenomes.withSRA.max0.5N.SNPs.forFbranch.100bp.$CHR.phased.vcf.gz) | gzip > allGenomes.withSRA.chr1-22.max0.5N.SNPs.forFbranch.100bp.$CHR.phased.headerFixed.vcf.gz
done

# Extract the individuals we want:
for CHR in chr{1..22}
do
bcftools view -S <(cut -f 1 LVRS.one.sample.per.group.info) \
    allGenomes.withSRA.chr1-22.max0.5N.SNPs.forFbranch.100bp.$CHR.phased.headerFixed.vcf.gz | \
    gzip > allGenomes.withSRA.max0.5N.SNPs.forFbranch.100bp.$CHR.phased.oneSamplePerGroup.vcf.gz
done

# Concatenate the chromosomes
vcf-concat allGenomes.withSRA.max0.5N.SNPs.forFbranch.100bp.chr{1..22}.phased.oneSamplePerGroup.vcf.gz | \
  gzip > allGenomes.withSRA.max0.5N.SNPs.forFbranch.100bp.chr1-22.phased.oneSamplePerGroup.vcf.gz

file=allGenomes.withSRA.max0.5N.SNPs.forFbranch.100bp.chr1-22.phased.oneSamplePerGroup


# Convert phased vcf to finestructure format
bsub < vcf2fineSTR.lsf

# Run finestructure
bsub -W 120:00 "fs ecoGroups.cp -phasefiles $prefix.chr{1..22}.SNPs.phased.fineSTRphase \
  -recombfiles $prefix.chr{1..22}.SNPs.phased.recomb -idfile fs.poplabels -go"
