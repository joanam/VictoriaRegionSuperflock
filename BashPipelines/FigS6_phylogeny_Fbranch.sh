
# Get a phylogeny including H. kimondo (piscivore from Lake Edward with very low sequencing depth)

module load gcc/4.8.2 gdc java/1.8.0_73 gatk/3.8-1

refgenome=/cluster/work/gdc/shared/p500/ref-genome/puncross.gapsEstimated.fasta
gatk=/cluster/apps/gatk/3.5/x86_64/GenomeAnalysisTK.jar

samples=""
for ind in `cut -f 1 edw_kivu.inds`
do
  samples=$samples" -V /cluster/work/gdc/shared/p500/victoriaGenomes/gvcf/"$ind.g.vcf.gz
done

for chr in chr{2..22}
do
  bsub -J gGVCF -R "rusage[mem=4000]" -W 24:00 -n 2  "java -Xmx12000M -Xms5000M -jar $gatk -T GenotypeGVCFs $samples -allSites -L $chr -R $refgenome -nt 2 -o $chr.edw.kivu.vcf.gz 2>&1 | tee $chr.edw.kivu.sterr.stdout.txt"
done

# Get the masks
cp /cluster/work/gdc/shared/p500/victoriaGenomes/masks/puncross.3masks.bed.gz .
gunzip puncross.3masks.bed.gz

# Filter the files and thin to 1 kb
module load gcc/4.8.2 gdc perl/5.18.4 zlib/1.2.8 vcftools/0.1.16

for chr in chr{1..22}
do
 bsub "vcftools --gzvcf edw_kivu_vic.$chr.vcf.gz --max-missing 0.5 --minDP 3 --thin 1000 --exclude-bed puncross.3masks.bed --stdout --recode | gzip > edw_kivu_vic.$chr.minDP3.max0.5N.masks.thin1kb.vcf.gz"
done

# Concatenate the vcf files
file="edw_kivu_vic.chr1-22.minDP3.max0.5N.masks.thin1kb"

bsub "vcf-concat edw_kivu_vic.chr{1..22}.minDP3.max0.5N.masks.thin1kb.vcf.gz | gzip > $file.vcf.gz"

# Then convert to phylip format for iqtree
module load python/2.7
bsub -J phylip "vcf2phylip.py -i $file.vcf.gz -o $file.phylip"

# Run IQtree (monomorphic sites included, so no ASC bias correction needed)
bsub -W 24:00 -R "rusage[mem=1500]" "iqtree2 -s $file.phylip -n 1000 -B 1000"


# Dsuite Fbranch

# First get a file of only SNPs
# Filter the files and thin to 1 kb
module load gcc/4.8.2 gdc perl/5.18.4 zlib/1.2.8 vcftools/0.1.16

for chr in chr{1..22}
do
 bsub "vcftools --gzvcf edw_kivu_vic.$chr.vcf.gz --remove-indels --max-alleles 2 --max-missing 0.5 --mac 2 --minDP 3 --thin 1000 --exclude-bed puncross.3masks.bed --stdout --recode | gzip > edw_kivu_vic.$chr.minDP3.max0.5N.masks.mac3.biallSNthin1kb.vcf.gz"
done

file="edw_kivu_vic.chr1-22.minDP3.max0.5N.masks.biallSNPs.thin1kb"
bsub "vcf-concat edw_kivu_vic.chr{1..22}.minDP3.max0.5N.masks.mac3.biallSNPs.thin1kb.vcf.gz | gzip > $file.vcf.gz"

# Run Dsuite on 100bp thinned dataset:
bsub -R "rusage[mem=7000]" -W 120:00 "Dsuite Dtrios -t fbranch.kimondo.tree.nwk -c -n kimondo \
$file.vcf.gz inds.fbranch.info.txt"

# Run Fbranch (Outgroup needs to be present in the tree)
module load gcc/4.9.2

bsub "Dsuite Dtrios -t fbranch.kimondo.tree.nwk -c -n kimondo edw_kivu_vic.chr1-22.minDP3.max0.5N.masks.thin1kb.vcf.gz inds.fbranch.info.txt"

Dsuite Fbranch -Z True fbranch.kimondo.tree.nwk inds.fbranch.info_kimondo_kimondo_tree.txt > fbranch.txt

Rscript `which removeNonsignFbranch.r` -i fbranch.txt -z 3 \
  -o fbranch.signOnly.txt

# Visualize it:
module load gdc python/3.6.1 zlib/1.2.8
dtools.py -n inds.fbranch.info.txt --use_distances --ladderize fbranch.signOnly.txt fbranch.kimondo.tree.nwk
