# Code, data and vectorised figures associated with the publication on the evolutionary origin of the Lake Victoria cichlid radiation.

## Initial data preparation
LSF job submission script to align the reads to the reference genome:
```sh
BashPipelines/00_alignToBam.lsf
```
LSF job submission script to use GATK Haplotype Caller to get gvcf files of each sample:
```sh
BashPipelines/00_runHaplotypeCaller.lsf
```
LSF job submission script to get a vcf file per chromosome and apply some basic filters
```sh
BashPipelines/00_runGenotypeGVCF_filter.sh
BashPipelines/00_removeTooHighDepthSites.sh # used to remove sites with too high sequencing depth, indicative of paralogous regions collapsed in the reference
```

## Analyses
Reconstruct and plot phylogenies 
```sh
# Generate the nuclear phylogeny
BashPipelines/01_iqtree2.sh
BashPipelines/01_vcf2phylip.py # used for file conversion

# Plot nuclear phylogenies as fan trees for Fig 2 and Figure S3
Rscripts/Fig2_FigS3_giantPhylogenies.R

# Generate the mitochondrial phylogeny
BashPipelines/01_iqtree2.sh

# Make co-phyloplots of quartet tree vs iqtree2
Rscripts/FigS15_cophyloplot_quartets_iqtree.R
```

Multidimensional Scaling (MDS) plots
```sh
Rscripts/Fig1_FigS4_FigS9_MDS.R
```

AMOVA: genomic variation explained by ecomorph or genus across the Victoria Radiation
```sh
BashPipelines/02_amova.sh
```

Admixture analyses
```sh
BashPipelines/Dstatistics.sh
BashPipelines/03_fbranch.sh # Example of Fbranch script, here for all LVRS and outgroups
Rscripts/Fig4B_f4values.R # plotting f4 values of ecomorph pairs across lakes
RScripts/FigS5_FigS8_Dstats.R # Plotting Dstats of Upper Nile ancestry in the LVRS and Western Lakes in the Victoria Radiation

```

