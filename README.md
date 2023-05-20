# Code, data and vectorised figures associated with the publication on the evolutionary origin of the Lake Victoria cichlid radiation


## Initial data preparation
* LSF job submission script to align the reads to the reference genome:
[BashPipelines/00_alignToBam.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/00_alignToBam.lsf)

* LSF job submission script to use GATK Haplotype Caller to get gvcf files of each sample:
[BashPipelines/00_runHaplotypeCaller.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/00_runHaplotypeCaller.lsf)

* LSF job submission script to get a vcf file per chromosome and apply some basic filters:
[BashPipelines/00_runGenotypeGVCF_filter.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/00_runGenotypeGVCF_filter.sh)

* Remove sites with too high sequencing depth, indicative of paralogous regions collapsed in the reference:
[BashPipelines/00_removeTooHighDepthSites.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/00_removeTooHighDepthSites.sh)


## Analyses

### Reconstruct and plot phylogenies

* Generate the nuclear phylogeny:
[BashPipelines/01_iqtree2.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/01_iqtree2.sh)

* for file conversion:
[BashPipelines/01_vcf2phylip.py](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/01_vcf2phylip.py)


* Plot nuclear phylogenies as fan trees for Fig 3 and Figure S2:
[Rscripts/Fig3_FigS2_giantPhylogenies.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/Fig3_FigS2_giantPhylogenies.R)


* Generate the mitochondrial phylogeny:
[BashPipelines/01_iqtree2.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/01_iqtree2.sh)


* Make co-phyloplots of quartet tree vs iqtree2:
[Rscripts/FigS14_cophyloplot_quartets_iqtree.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/FigS14_cophyloplot_quartets_iqtree.R)


* Co-phyloplot of nuclear vs mitochondrial phylogenies for Aux Fig. 1
[Rscripts/AuxFig1_cophyloplot_mtDNA-nuc.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/AuxFig1_cophyloplot_mtDNA-nuc.R)


### Multidimensional Scaling (MDS) plots
Compute MDS with different combinations of individuals: [Rscripts/FigS3_FigS8_MDS.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/FigS3_FigS8_MDS.R)



### AMOVA
Compute genomic variation explained by ecomorph or genus across the Victoria Radiation:
[BashPipelines/02_amova.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/02_amova.sh)



### Admixture analyses
* Compute D statistics with ADMIXTOOLS: [BashPipelines/Dstatistics.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/Dstatistics.sh)

* Example of Fbranch script, here for all LVRS and outgroups:
[BashPipelines/03_fbranch.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/03_fbranch.sh)

* Plotting f4 values of ecomorph pairs across lakes:
[Rscripts/Fig5B_f4values.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/Fig5B_f4values.R)

* Plotting Dstats of Upper Nile ancestry in the LVRS and Western Lakes in the Victoria Radiation:
[RScripts/FigS4_FigS7_Dstats.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/FigS4_FigS7_Dstats.R)

* Compute DFOIL to test the direction of gene flow between Har. squamipinnis (Edward piscivore) and Victoria piscivores:
[BashPipelines/FigS6_computeDfoil.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS6_computeDfoil.sh)
which requires
[BashPipelines/FigS6_runDfoil.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS6_runDfoil.sh)

* Compute fdM to test for correlations between ancestries (example of dwarf predators): [FigS13_fdM_dwarfPredators.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS13_fdM_dwarfPredators.sh)


### Ecomorph associated alleles

* Get allele frequencies for each ecomorph and group:
[BashPipelines/Fig5_FigS10-S12_getFrqs.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/Fig5_FigS10-S12_getFrqs.lsf)

* Get sites with at least 0.9 allele frequency difference between ecomorph and all other Victoria Radiation cichlids:
[BashPipelines/Fig5_FigS10-S12_getDiffSites.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/Fig5_FigS10-S12_getDiffSites.lsf)

* Plot them:
[RScripts/FigS10-S12_alleleFrequencies_ecogroups.R](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/FigS10-S12_alleleFrequencies_ecogroups.R)



### DXY between allopatric and sympatric piscivores and paedophages
* Get the file with all sites of piscivores and paedophages (one sample per species):
[FigS9_runGenotypeGVCF.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS9_runGenotypeGVCF.lsf)

* Compute DXY:
[FigS9_computeDXY.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS9_computeDXY.sh)

* Plot DXY:
[RScripts/FigS9_dxy_sympatric_vs_allopatric.pdf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/FigS9_dxy_sympatric_vs_allopatric.pdf)



### FineSTRUCTURE

Run fineSTRUCTURE: [BashPipelines/FigS6_runFineStructure.sh](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS6_runFineStructure.sh)

Requires these scripts:

* [BashPipelines/FigS6_vcf2fineSTR.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS6_vcf2fineSTR.lsf)
* [BashPipelines/FigS6_runBeagle.lsf](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/BashPipelines/FigS6_runBeagle.lsf)
* [RScripts/FigS6_addRecombRates.r](https://github.com/joanam/VictoriaRegionSuperflock/blob/main/Rscripts/FigS6_addRecombRates.R)


## Additional files:

* Figures in vector format can be found [here](https://github.com/joanam/VictoriaRegionSuperflock/tree/main/Figures)

* The phylogenies can be found [here](https://github.com/joanam/VictoriaRegionSuperflock/tree/main/Phylogenies)
