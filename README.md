# Code, data and vectorised figures associated with the publication on the evolutionary origin of the Lake Victoria cichlid radiation


## Initial data preparation
### LSF job submission script to align the reads to the reference genome:
[BashPipelines/00_alignToBam.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/00_alignToBam.lsf)

### LSF job submission script to use GATK Haplotype Caller to get gvcf files of each sample:
[BashPipelines/00_runHaplotypeCaller.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/00_runHaplotypeCaller.lsf)

### LSF job submission script to get a vcf file per chromosome and apply some basic filters:
[BashPipelines/00_runGenotypeGVCF_filter.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/00_runGenotypeGVCF_filter.sh)

Remove sites with too high sequencing depth, indicative of paralogous regions collapsed in the reference:
[BashPipelines/00_removeTooHighDepthSites.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/00_removeTooHighDepthSites.sh)


## Analyses

### Reconstruct and plot phylogenies 

* Generate the nuclear phylogeny:
[BashPipelines/01_iqtree2.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/01_iqtree2.sh)

* for file conversion:
[BashPipelines/01_vcf2phylip.py](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/01_vcf2phylip.py)


* Plot nuclear phylogenies as fan trees for Fig 2 and Figure S3:
[Rscripts/Fig2_FigS3_giantPhylogenies.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/Fig2_FigS3_giantPhylogenies.R)


* Generate the mitochondrial phylogeny:
[BashPipelines/01_iqtree2.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/01_iqtree2.sh)


* Make co-phyloplots of quartet tree vs iqtree2:
[Rscripts/FigS15_cophyloplot_quartets_iqtree.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/FigS15_cophyloplot_quartets_iqtree.R)



### Multidimensional Scaling (MDS) plots
Compute MDS with different combinations of individuals: [Rscripts/Fig1_FigS4_FigS9_MDS.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/Fig1_FigS4_FigS9_MDS.R)



### AMOVA
Compute genomic variation explained by ecomorph or genus across the Victoria Radiation: 
[BashPipelines/02_amova.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/02_amova.sh)



### Admixture analyses
* Compute D statistics with ADMIXTOOLS: [BashPipelines/Dstatistics.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/Dstatistics.sh)

* Example of Fbranch script, here for all LVRS and outgroups: 
[BashPipelines/03_fbranch.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/03_fbranch.sh)

* Plotting f4 values of ecomorph pairs across lakes: 
[Rscripts/Fig4B_f4values.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/Fig4B_f4values.R)

* Plotting Dstats of Upper Nile ancestry in the LVRS and Western Lakes in the Victoria Radiation: 
[RScripts/FigS5_FigS8_Dstats.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/FigS5_FigS8_Dstats.R)

* Compute DFOIL to test the direction of gene flow between Har. squamipinnis (Edward piscivore) and Victoria piscivores: 
[BashPipelines/FigS7_computeDfoil.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS7_computeDfoil.sh)
which requires 
[BashPipelines/FigS7_runDfoil.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS7_runDfoil.sh)

* Compute fdM to test for correlations between ancestries (example of dwarf predators): [FigS14_fdM_dwarfPredators.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS14_fdM_dwarfPredators.sh)


### Ecomorph associated alleles

* Get allele frequencies for each ecomorph and group:
[BashPipelines/Fig4_FigS11-S13_getFrqs.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/Fig4_FigS11-S13_getFrqs.lsf)

* Get sites with at least 0.9 allele frequency difference between ecomorph and all other Victoria Radiation cichlids: 
[BashPipelines/Fig4_FigS11-S13_getDiffSites.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/Fig4_FigS11-S13_getDiffSites.lsf)

* Plot them: 
[RScripts/FigS11-S13_alleleFrequencies_ecogroups.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/S13_alleleFrequencies_ecogroups.R)



### DXY between allopatric and sympatric piscivores and paedophages
* Get the file with all sites of piscivores and paedophages (one sample per species):
[FigS10_runGenotypeGVCF.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS10_runGenotypeGVCF.lsf)

* Compute DXY: 
[FigS10_computeDXY.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS10_computeDXY.sh)

* Plot DXY: 
[RScripts/FigS10_dxy_sympatric_vs_allopatric.pdf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/FigS10_dxy_sympatric_vs_allopatric.pdf)



### FineSTRUCTURE

Run fineSTRUCTURE: [BashPipelines/FigS7_runFineStructure.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS7_runFineStructure.sh)

Requires these scripts:

* [BashPipelines/FigS7_vcf2fineSTR.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS7_vcf2fineSTR.lsf)
* [BashPipelines/FigS7_runBeagle.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS7_runBeagle.lsf)
* [RScripts/FigS7_addRecombRates.r](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/FigS7_addRecombRates.R)


## Additional files:

* Figures in vector format can be found [here](https://github.com/joanam/VictoriaRegionSuperflock/tree/main/Figures)

* The phylogenies can be found [here](https://github.com/joanam/VictoriaRegionSuperflock/tree/main/Phylogenies)

