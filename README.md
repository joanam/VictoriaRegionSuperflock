# Code, data and vectorised figures associated with the publication on the evolutionary origin of the Lake Victoria cichlid radiation.

# Initial data preparation
LSF job submission script to align the reads to the reference genome:
```sh
00_alignToBam.lsf
```
LSF job submission script to use GATK Haplotype Caller to get gvcf files of each sample:
```sh
00_runHaplotypeCaller.lsf
```
LSF job submission script to get a vcf file per chromosome and apply some basic filters
```sh
00_runGenotypeGVCF_filter.sh
00_removeTooHighDepthSites.sh # used to remove sites with too high sequencing depth, indicative of paralogous regions collapsed in the reference
```

# Generate phylogenies with iqtree2
```sh
01_iqtree2.sh
01_vcf2phylip.py # used for file conversion
```


