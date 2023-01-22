# Dsuite with species of Victoria radiation all separately:

file=/cluster/work/gdc/shared/p500/victoriaGenomes/LVRS/newLVRS/vcf/allGenomes.withSRA.chr1-22.max0.5N.SNPs.vcf.gz

# Filter the dataset by extracting the subset of individuals of interest
# Prune to 1 SNP per 500 bp
module load gcc/4.8.2 gdc perl/5.18.4 zlib/1.2.8 vcftools/0.1.16
bsub -R "rusage[mem=3000]" -W 120:00 "vcftools --gzvcf $file \
 --keep species.forFbrach.info --thin 500 --mac 2 \
 --max-missing 0.75 --recode --stdout | gzip > allGenomes.withSRA.chr1-22.max0.5N.SNPs.forFbranch.500bp.vcf.gz"

# Get sample with least missing data of each group
file=allGenomes.withSRA.chr1-22.max0.5N.SNPs.forFbranch.500bp.vcf.gz
bsub -R "rusage[mem=3000]" -W 120:00 "vcftools --gzvcf $file \
 --keep species.withLVRS.forFbrach.info --missing-indv \
 --out species.withLVRS.forFbrach."

paste <(awk 'NR==FNR{a[$1]=$2; next}{print $0"\t"a[$1]}' species.withLVRS.forFbrach.info allGenomes.withSRA.chr21.imiss) | \
 sort -nrk 5 | awk '{a[$6]=$1"\t"$6}END{for(i in a) print a[i]}' > one.sample.per.group.txt


# Extract sample with least missing data of each group to make a tree
module load gcc/4.8.2 gdc perl/5.18.4 zlib/1.2.8 vcftools/0.1.16

bsub -R "rusage[mem=3000]" -W 120:00 "vcftools --gzvcf $file \
 --keep one.sample.per.group.txt --recode --stdout | \
   gzip > one.sample.per.group.chr1-21.thin500bp.vcf.gz"

# Convert to phylip and run iqtree2 to get a backbone tree for fbranch
file=one.sample.per.group.chr1-21.thin500bp

bsub "vcf2phylip.py -i $file.vcf.gz -o $file.phylip"

iqtree2 -s $file.phylip -m GTR+ASC -T 1 -o 131282 --prefix $file
bsub -W 120:00 -R "rusage[mem=1000]" -n 16 "iqtree2 -redo -s $file.varsites.phy -st DNA -m GTR+ASC -T 16 -o 131282 --prefix LVRS.one.sample.per.group"

# Generate a script that replaces the tree labels by group name and run it
awk -v q="'" 'BEGIN{FS="\t"; print "#!/bin/bash\ncat $1 | sed \\"}
{print "-e "q"s/"$1":/"$1"_"$2":/g"q" \\""-e "q"s/"$1")/"$1"_"$2")/g"q" \\"}
END{print " > $1.renamed"}' species.withLVRS.forFbrach.info > rename.sh
chmod a+x rename.sh # make it executable
rename.sh LVRS.one.sample.per.group.phylip.varsites.phy.treefile

# Run Dsuite Dtrios to get Dstats and f4 values
# Note, the files are renamed to be a bit more informative

tree="species.tree.for.Fbranch.burtOut.nwk"
sampleInfo="species.forFbranch.burtOut.info"
prefix="lvrs.speciesTree.burtOut"
vcf="allGenomes.withSRA.chr1-22.max0.5N.SNPs.forFbranch.500bp.vcf.gz"

bsub -R "rusage[mem=40000]" -W 120:00 "Dsuite Dtrios -t $tree \
 -c -n $prefix $vcf $sampleInfo"


# Run fbranch with zscores
module load gcc/4.9.2

bsub "Dsuite Fbranch -Z True $tree $prefix_tree.txt > fbranch.Zscores.txt"


# Adjust non-significant D stats to 0
Rscript `which removeNonsignFbranch.r` -i fbranch.Zscores.txt -z 3 \
  -o fbranch.signOnly.txt


 # Visualize it the fbranch result
 module load gdc python/3.6.1 zlib/1.2.8
 dtools.py -n $prefix --use_distances --color-cutoff 1 --tree-label-size 18 \
  --ladderize fbranch.signOnly.txt $tree
