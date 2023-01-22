path=/cluster/home/meierjo/bin/

module load gatk/3.5

for i in {1..22}
do
bsub -J gGVCF -R "rusage[mem=4000]" -q "normal.120h" -n 4 \
 "java -Xmx16000M -Xms5000M \
  -jar /cluster/apps/gatk/3.5/x86_64/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs $samples -allSites \
    -R $ref -nt 4 -L chr$i \
    -o /cluster/scratch/meierjo/vcf/piscivores_paedophages_VictoriaKyoga.chr$i.vcf.gz \
    2>&1 | tee GenotypeVCF.stdout.stderr.txt"
done

# Run pixy to compute dxy across all pairs of piscivores and paedophages
pixy --stats dxy  --vcf piscivores_paedophages_VictoriaKyoga.chr1.max0.5N.minDP5.vcf.gz  --populations populations.info  --window_size 1000000000  --n_cores 2

