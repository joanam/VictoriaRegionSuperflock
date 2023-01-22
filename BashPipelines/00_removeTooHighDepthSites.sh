#! /bin/bash
# Usage: removeTooHighDepthSites.sh <file.vcf> <maxDepth>
# Creates a file with all sites failing the maxDepth filter and a new vcf file $1.max$2meanSiteDepth
# Created by Joana Meier, Dec-21, 2015

# This clause checks if there were two arguments given
if [ $# -lt 2 ]
then
        echo -e "ERROR: Not enough arguments provided!\nUsage: removeTooHighDepthSites.sh <file.vcf> <maxDepth>"
        exit 1
fi

# This clause checks if the VCF file exists
if [ -f $1 ]
then
        :
else
        echo -e "ERROR: VCF file does not exist!\nUsage: removeTooHighDepthSites.sh <file.vcf> <maxDepth>>"
        exit 1
fi

# Get bed file with sites failing the max-depth filter
vcftools15 --vcf $1 --site-mean-depth --stdout | \
 awk -v maxD=$2 'NR>1 {if($3>maxD) print $1"\t"$2}' > ${1%.vcf}.tooHighDepthSites

if [ -s ${1%.vcf}.tooHighDepthSites ]
then
	# remove the sites with too High Depth:
	vcftools15 --vcf $1 --exclude-positions ${1%.vcf}.tooHighDepthSites \
		--recode --out ${1%.vcf}.max$2meanSiteDepth
	repairmissing.sh ${1%.vcf}.max${2}meanSiteDepth.recode.vcf
	echo "Finished! Filtered out "`cat ${1%.vcf}.tooHighDepthSites | wc -l`" sites"
	echo "New file without these sites: ${1%.vcf}.max${2}meanSiteDepth.vcf"
else
	echo -e "ERROR: No sites failed the maxDepth filter."
	rm ${1%.vcf}.tooHighDepthSites
fi
