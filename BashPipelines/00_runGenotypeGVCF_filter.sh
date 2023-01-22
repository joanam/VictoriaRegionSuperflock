# Get vcf of chr1
module load gcc/4.8.2 gdc java/1.8.0_73 gatk/3.8-1

refgenome=/cluster/work/gdc/shared/p500/ref-genome/puncross.gapsEstimated.fasta
gatk=/cluster/apps/gatk/3.5/x86_64/GenomeAnalysisTK.jar

# Get all the samples with preceeding -V into a string variable
samples="";
for ind in `cut -f 1 edw_kivu.inds`;
do
   samples=$samples" -V /cluster/work/gdc/shared/p500/victoriaGenomes/gvcf/"$ind.g.vcf.g
done

# Submit each chromosome as separate job:
for chr in chr{2..22}
do
bsub -J gGVCF -R "rusage[mem=4000]" -W 24:00 -n 2  "java -Xmx12000M -Xms5000M -jar $gatk -T GenotypeGVCFs $samples \
   -allSites -L $chr -R $refgenome -nt 2 -o $chr.vcf.gz 2>&1 | tee $chr.sterr.stdout.txt"
done

# Get the masks
cp /cluster/work/gdc/shared/p500/victoriaGenomes/masks/puncross.3masks.bed.gz .
gunzip puncross.3masks.bed.gz

# Filter the files and thin to 1 kb
module load gcc/4.8.2 gdc perl/5.18.4 zlib/1.2.8 vcftools/0.1.16

for chr in chr{1..22}
do
 bsub "vcftools --gzvcf $chr.vcf.gz --max-missing 0.5 --minDP 6 --minGQ 20 --remove-indels \
  --exclude-bed puncross.3masks.bed --stdout --recode | gzip > $chr.filtered.vcf.gz"
done

# Remove sites with too high depth:
for chr in chr{1..22}
do
  bsub "removeTooHighDepthSites.sh  $chr.filtered.vcf.gz <maxDepth>"
done
