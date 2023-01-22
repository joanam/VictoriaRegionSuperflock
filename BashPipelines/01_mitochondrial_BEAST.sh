#!/bin/bash

# Mitochondrial BEAST phylogeny

# Download samples from Svardal et al. (potential Eastern clade taxa and two Malawi reps)

module load new gcc/5.2.0 gdc sra-tools/2.10.8
fastq-dump --split-files --gzip ERR715501   # SAMEA2661249 Astatotilapia sp. 'Ruaha' from Rujewa
fastq-dump --split-files --gzip ERR1822251  # SAMEA4033334 Astatotilapia bloyeti from Lake Kumba
fastq-dump --split-files --gzip ERR1822250  # SAMEA4033333 Astatotilapia bloyeti from Lake Burungi
fastq-dump --split-files --gzip ERR3634113  # SAMEA6147772 A. Ruaha "red cheek" from Kidatu, Gt Ruaha, Rufiji System
fastq-dump --split-files --gzip ERR1081365  # SAMEA3388853 Rhamphochromis longiceps
fastq-dump --split-files --gzip ERR1081384  # SAMEA3388874 Maylandia zebra
fastq-dump --split-files --gzip ERR715538   # SAMEA2661389 A. calliptera Near_Kyela Malawi catchment

# Move the files to the fastq group folder renaming them to fit naming standards (.1, not _1)
for i in *_2.fastq.gz
do
  cp $i /cluster/work/gdc/shared/p500/victoriaGenomes/fastq/${i%_2.fastq.gz}.2.fastq.gz
  cp ${i%_2.fastq.gz}_1.fastq.gz /cluster/work/gdc/shared/p500/victoriaGenomes/fastq/${i%_2.fastq.gz}.1.fastq.gz
done

# Align them to our P. nyererei genome:
ls *_2.fastq.gz | sed 's/_2.fastq.gz//' > individuals
bsub < `which alignToBam.lsf`

# Get gvcf files:
bsub < `which runHaplotypeCaller.lsf`

# Get a vcf file of only the mitochondrial genome
module load gcc/4.8.2 gdc java/1.8.0_73 gatk/3.5

# Get a list of the gvcf files to include
path=/cluster/work/gdc/shared/p500/victoriaGenomes/gvcf/

samples=""
for i in `cut -f 1 samples.forAstral.Balanced.txt`
do
 samples=$samples" -V "$path$i.g.vcf.gz
done

# add new samples from SRA:
for ind in `cat individuals.SRA.txt`
do
  samples=$samples" -V "$ind.g.vcf.gz
done

# Load the module required for GATK
module load gcc/4.8.2 gdc java/1.8.0_73 gatk/3.8-1

# Specify the reference genome and path to GATK jar file
refgenome=/cluster/work/gdc/shared/p500/ref-genome/puncross.gapsEstimated.fasta
gatk=/cluster/apps/gatk/3.5/x86_64/GenomeAnalysisTK.jar

# Run GenotypeGVCF
bsub -J gGVCF -R "rusage[mem=3000]" -W 24:00 -n 4 \
 "java -Xmx12000M -Xms5000M \
  -jar $gatk -T GenotypeGVCFs $samples -allSites -L chrM \
    -R $refgenome -nt 4 -o mtDNA.AstralSamples.Balanced.vcf.gz \
    2>&1 | tee mtDNA.GenotypeVCF.stdout.stderr.txt"

file=mtDNA.AstralSamples.Balanced

# Some filtering:
module load gcc/4.8.2 gdc perl/5.18.4 zlib/1.2.8 vcftools/0.1.15
bsub "vcftools --gzvcf $file.vcf.gz --minDP 10 \
  --recode --stdout | gzip > $file.minDP10.vcf.gz"

# Convert vcf file to beast2 input xml file
file=mtDNA.AstralSamples.Balanced.minDP10
vcf2phylip.py -i $file.vcf.gz -o $file.phylip -f # -f fills missing sites with Ns
awk 'NR>1 {print ">"$1"\n"$2}' $file.phylip  > $file.fasta

# Run partition finder
# Manually make: partition_finder.cfg # Following the manual

module load gcc/4.8.2 gdc hdf5/1.8.12 open_mpi/1.6.5 boost/1.63.0 python/2.7.6
bsub -W 120:00 "python partitionfinder-2.1.1/PartitionFinder.py --no-ml-tree PartitionFinder/"


# Download and load into Beauti

# Modify file by hand to create a nexus file showing the partitions as from PartitionFinder

# In Beauti: import alignment, select all of partitions (text, no clicking of boxes) and link trees then clocks
# Click multiple times on link trees and clocks until they all show the same value
# Rename them to tree and clock
# For each partition select site model as BEAST model test
# Choose Calibrated Yule Model
# set a normal prior on the Malawi cichlid crown and the modern haplochromine stem

# Note, I downloaded and unpacked BEAST and then I had to add BEASTLabs and bModelTest with the packagemanager
module load gcc/4.8.2 gdc java/1.8.0_101 beast2/2.4.6 beagle/4.1

# Run BEAST with different calibration sets

bsub -W 24:00 -n 6 -R "rusage[mem=400]" \
 "beast/bin/beast -beagle_CPU -threads 8 Astral.samples.balanced.cichlidFossils.alluaOutgroup.xml"

bsub -W 24:00 -n 6 -R "rusage[mem=400]" \
  "beast/bin/beast -beagle_CPU -threads 8 Astral.samples.balanced.nonCichlidFossils.alluaOut.xml"

bsub -W 24:00 -n 6 -R "rusage[mem=400]" \
  "beast/bin/beast -beagle_CPU -threads 8 Astral.samples.balanced.C10.alluaOut.xml"

bsub -W 24:00 -n 6 -R "rusage[mem=300]" \
  "beast/bin/beast -beagle_CPU -threads 8 Astral.samples.balanced.Gondwana.alluaOut.xml"

bsub -W 24:00 -n 6 -R "rusage[mem=400]" \
    "beast/bin/beast -beagle_CPU -threads 8 Astral.samples.balanced.lacPaleo.alluaOut.xml"
