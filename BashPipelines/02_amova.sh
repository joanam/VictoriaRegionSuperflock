# AMOVA with poppr

# Extract samples of species with at least three individuals for AMOVA
module load gcc/4.8.2 gdc perl/5.18.4 zlib/1.2.8 vcftools/0.1.16
for chr in chr{1..22}; do bsub -R "rusage[mem=3000]" -J "extr"$chr "vcftools --gzvcf allGenomes.withSRA.$chr.max0.5N.SNPs.vcf.gz --keep samples.for.AMOVA --max-missing 0.75 --minDP 5 --recode --mac 1 --stdout | gzip > allGenomes.withSRA.$chr.max0.25N.minDP5.AMOVA.samples.vcf.gz"; done
for chr in chr{1..22}; do bsub -R "rusage[mem=3000]" -J "prune"$chr  -w "done(extrchr2)" "ldPruning.sh allGenomes.withSRA.$chr.max0.25N.minDP5.AMOVA.samples.vcf.gz 0.05"; done
bsub -R "rusage[mem=5000]" -w "done(prunechr2)" "vcf-concat allGenomes.withSRA.chr{1..22}.max0.25N.minDP5.AMOVA.samples.LDpruned.vcf.gz | gzip > allGenomes.withSRA.chr1-22.max0.25N.minDP5.AMOVA.samples.LDpruned.vcf.gz"


# Again for all victoria-kyoga radiation samples (no satellite lakes)
module load gcc/4.8.2 gdc perl/5.18.4 zlib/1.2.8 vcftools/0.1.16
path=/cluster/work/gdc/shared/p500/victoriaGenomes/LVRS/newLVRS/vcf
for chr in chr{1..22}; do bsub -R "rusage[mem=3000]" -J "extr"$chr "vcftools --gzvcf ${path}/allGenomes.withSRA.$chr.vcf.gz --keep victoria.samples --max-missing 0.75 --minDP 5 --recode --mac 2 --stdout | gzip > allGenomes.withSRA.$chr.max0.25N.minDP5.VicKyo.samples.vcf.gz"; done
for chr in chr{1..22}; do bsub -R "rusage[mem=3000]" -J "prune"$chr  -w "done(extrchr2)" "ldPruning.sh allGenomes.withSRA.$chr.max0.25N.minDP5.VicKyo.samples.vcf.gz 0.05"; done
bsub -R "rusage[mem=5000]" -w "done(prunechr2)" "vcf-concat allGenomes.withSRA.chr{1..22}.max0.25N.minDP5.VicKyo.samples.LDpruned.vcf.gz | gzip > allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.LDpruned.vcf.gz"


# Check that none of the individuals has high missing data proportion:
vcftools --gzvcf allGenomes.withSRA.chr1-22.max0.25N.minDP5.AMOVA.samples.LDpruned.vcf.gz --missing-indv
vcftools --gzvcf allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.LDpruned.vcf.gz --missing-indv

# Get site depth information
vcftools --gzvcf allGenomes.withSRA.chr1-22.max0.25N.minDP5.AMOVA.samples.LDpruned.vcf.gz --site-mean-depth
vcftools --gzvcf allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.LDpruned.vcf.gz --site-mean-depth


# Get info file in the right order
awk 'FNR==NR {x2[$1] = $0; next} $1 in x2 {print x2[$1]}' samples.for.AMOVA.info out.imiss > samples.for.AMOVA.info.tmp
mv samples.for.AMOVA.info.tmp samples.for.AMOVA.info

awk 'FNR==NR {x2[$1] = $0; next} $1 in x2 {print x2[$1]}' \
 victoria.info allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.LDpruned.imiss \
  > samples.for.AMOVA.VicKyo.info

# Remove sites with particularly high or low sequencing depth
R
depths<-read.table("allGenomes.withSRA.chr1-22.max0.25N.minDP5.AMOVA.samples.LDpruned.ldepth.mean",header=T)
depths<-read.table("allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.LDpruned.ldepth.mean",header=T)

summary(depths$MEAN_DEPTH)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   6.375   26.560   28.420   28.208   29.710 4294.260
write.table(file="allGenomes.withSRA.chr1-22.max0.25N.minDP5.AMOVA.samples.depth20-40.sites",
  x=depths[depths$MEAN_DEPTH>20 & depths$MEAN_DEPTH<40,1:2],quote=F,row.names=F)
write.table(file="allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.depth20-40.sites",
    x=depths[depths$MEAN_DEPTH>20 & depths$MEAN_DEPTH<40,1:2],quote=F,row.names=F)
quit()

bsub "vcftools --gzvcf allGenomes.withSRA.chr1-22.max0.25N.minDP5.AMOVA.samples.LDpruned.vcf.gz \
  --positions allGenomes.withSRA.chr1-22.max0.25N.minDP5.AMOVA.samples.depth20-40.sites \
 --recode --stdout | gzip > allGenomes.withSRA.chr1-22.max0.25N.minDP5.AMOVA.samples.LDpruned.depth20-40.vcf.gz"

bsub "vcftools --gzvcf allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.LDpruned.vcf.gz \
   --positions allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.depth20-40.sites \
  --recode --stdout | gzip > allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.LDpruned.depth20-40.vcf.gz"

# file with only n>=3 species: We remain with 3,906,147 sites out of 4,248,008 sites
# file with all VicKyo species: 2,232,562 out of 2,509,595 sites remain

# Randomly subsample 1 mio sites
module load gcc/4.8.2 gdc vcflib/1.0.1.1
bsub -R "rusage[mem=3000]" "vcfrandomsample -r 0.256 \
 <(zcat allGenomes.withSRA.chr1-22.max0.25N.minDP5.AMOVA.samples.LDpruned.depth20-40.vcf.gz) | \
  gzip > allGenomes.withSRA.chr1-22.max0.25N.minDP5.AMOVA.samples.LDpruned.depth20-40.1mioSites.vcf.gz"

module load gcc/4.8.2 gdc vcflib/1.0.1.1
bsub -R "rusage[mem=3000]" "vcfrandomsample -r 0.45 \
   <(zcat allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.LDpruned.depth20-40.mac2.noNA.vcf.gz) | \
    gzip > allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.LDpruned.depth20-40.1mioSites.vcf.gz"


# Compute AMOVA in an interactive node
bsub -R "rusage[mem=10000]" -W 24:00 -Is /bin/bash
module load new gcc/6.3.0 gdc udunits zlib/1.2.8 r/3.6.0

R

# Load the required packages
require('poppr')
library('adegenet')
library('vcfR')

# Read in the vcf file
#vcf <- read.vcfR('allGenomes.withSRA.chr1-22.max0.25N.minDP5.AMOVA.samples.LDpruned.depth20-40.1mioSites.vcf.gz')
vcf <- read.vcfR('allGenomes.withSRA.chr1-22.max0.25N.minDP5.VicKyo.samples.LDpruned.depth20-40.1mioSites.vcf.gz')
geno <- vcfR2genlight(vcf)

# Add species and genus information
#info<-read.table('samples.for.AMOVA.info',header=T)
info<-read.table('samples.for.AMOVA.VicKyo.info',header=T)
strata(geno)<-info

# Compute AMOVA
amovaGenusSp<-poppr.amova(geno, ~Genus/Species)
amovaEcoSp<-poppr.amova(geno, ~Diet/Species)

# Results:
write.table(file="amovaGenusSp.allSpecies.variance",amovaGenusSp$componentsofcovariance)
write.table(file="amovaGenusSp.allSpecies.results",amovaGenusSp$results)

write.table(file="amovaEcoSp.allSpecies.variance",amovaEcoSp$componentsofcovariance)
write.table(file="amovaEcoSp.allSpecies.results",amovaEcoSp$results)

                            "Sigma" "%"
"Variations  Between Diet  " 627.543896801162 2.49745567791712
"Variations  Between Species Within Diet" 1815.98491849276 7.22713083310987
"Variations  Between samples Within Species" 70.4559785037618 0.280395816857382
"Variations  Within samples" 22613.3438851454 89.9950176721156

"Variations  Between Genus  " 803.894560387341 3.19241536306459
"Variations  Between Species Within Genus" 1677.11078448712 6.66011999307254
"Variations  Between samples Within Species" 87.0406883531959 0.345654821419085
"Variations  Within samples" 22613.3438851454 89.8018098224438

# Genus explains 3.192 %, species within genus: 6.660%, whereas diet only explains 2.497%

# Test if the species are significantly different
set.seed(1999)
GenusSpSign   <- randtest(amovaGenusSp, nrepet = 999)
EcoSpSign   <- randtest(amovaEcoSp, nrepet = 999)

pdf(file="randomTestANOVA.pdf")
plot(GenusSpSign)
plot(EcoSpSign)
dev.off()

quit()

# Weirdly the plot is a pdf, called R plot. Renamed:
cp Rplots.pdf GenusEcologySpecies.comparisonToRandom.pdf

# p-values of variation between species, diet and genus all highly significant p<0.01
