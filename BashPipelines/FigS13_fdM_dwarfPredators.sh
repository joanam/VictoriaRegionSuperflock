
# Compute fdM to study the admixed ancestry of dwarf predators


 # Extract the dwarf allies from the big vcf files:
 for chr in chr{1..22}
 do
  bsub -J dwarf$chr -R "rusage[mem=3000]" -W 120:00  "vcftools --gzvcf allGenomes.withSRA.$chr.vcf.gz \
     --max-missing 0.5 --mac 1 --remove-indels --recode --stdout --keep dwarf.inds | \
     gzip > dwarfAllies.$chr.vcf.gz"
 done

 bsub "vcf-concat dwarfAllies.chr{1..22}.vcf.gz | gzip > dwarfAllies.chr1-22.vcf.gz"

 # Some more stringent filtering
 bsub -J filter "vcftools --gzvcf dwarfAllies.chr1-22.vcf.gz --minDP 5 --max-missing 0.75 --recode --stdout | \
  gzip > dwarfAllies.chr1-22.minDP5.max0.25N.vcf.gz"

 # LD pruning
 bsub -J prune -w "done(filter)" "ldPruning.sh dwarfAllies.chr1-22.minDP5.max0.25N.vcf.gz"


 # fd zoopl into dwarfPred
 module load python
 for i in $file.chr{1..22}
 do
  bsub  -J "fd"$i -W 24:00 -n 1 \
  "python ~/bin/ABBABABAwindows.py \
    -w 50000 -m 10 -s 50000 -g $file.geno.gz  \
    -o $file.50kb.pisc-dwarf-zoopl-allua.csv -f phased -T 1 \
    --minData 0.5 --writeFailedWindows \
    -P1 pisc 130704,130732,130712,Ig279,checkmate,11050,13135,161543,130735,130705 \
    -P2 dwarfP 109846,13819,Ma41,Ma71,104621,104689,79547,13926,79542,103778 \
    -P3 zoopl 79559,79560,79561,IG104,103754,Y_P,109761 \
    -O out 105845,103218,131201"
 done

 # Combine the results of the different chromosomes
 head -1 $file.chr1.50kb.pisc-dwarf-zoopl-allua.csv > $file.50kb.pisc-dwarf-zoopl-allua.csv
 cat $file.chr{1..22}.50kb.pisc-dwarf-zoopl-allua.csv | \
  grep -v scaffold >> $file.50kb.pisc-dwarf-zoopl-allua.csv


 # fd pisc into dwarfPred
 for file in $file.chr{1..22}
 do
  bsub  -J "fd"$i -W 24:00 -n 1 \
  "python ~/bin/ABBABABAwindows.py \
    -w 50000 -m 10 -s 50000 -g $file.geno.gz  \
    -o $file.50kb.zoopl-dwarf-pisc-allua.csv -f phased -T 1 \
    --minData 0.5 --writeFailedWindows \
    -P3 pisc 130704,130732,130712,Ig279,checkmate,11050,13135,161543,130735,130705 \
    -P2 dwarfP 109846,13819,Ma41,Ma71,104621,104689,79547,13926,79542,103778 \
    -P1 zoopl 79559,79560,79561,IG104,103754,Y_P,109761 \
    -O out 105845,103218,131201"
 done


 # Combine the results of the different chromosomes
 head -1 $file.chr1.50kb.zoopl-dwarf-pisc-allua.csv > $file.50kb.zoopl-dwarf-pisc-allua.csv
 cat $file.chr{1..22}.50kb.zoopl-dwarf-pisc-allua.csv | \
  grep -v scaffold >> $file.50kb.zoopl-dwarf-pisc-allua.csv
