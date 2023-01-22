#!/bin/bash

# Get a phylogeny with iqtree2

# Convert vcf to phylip
bsub "vcf2phylip.py -i allGenomes.1kbThinned.lakes.tree.groups.withSRA.chr1-21.noSquami.vcf.gz -o allGenomes.1kbThinned.lakes.tree.groups.withSRA.chr1-21.noSquami.phylip"

# Example file
file=allGenomes.1kbThinned.withSRA.chr1-21.phylip

# Get phylip file with sites that are not polymorphic (ascertainment bias correction as only SNPs in this file)
iqtree2 -s $file -m MFP+ASC

# Run IQtree2 with 1000 rapid bootstraps
bsub -n 8 -W 120:00 "iqtree2 -s $file.varsites.phy -B 1000 -bnni -m MFP+ASC -nt 8"
