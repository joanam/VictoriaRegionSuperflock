# Run Dfoil to infer the direction of gene flow with Har. squamipinnis (Edward predator)

module load gcc/4.8.2 gdc perl/5.18.4 zlib/1.2.8 vcftools python/2.7; perl-init

prefix=nubila-VicPisc-squami-KivPisc-burtoni

# run script runDfoil.sh for all combinations of Har. squamipinnis, Victoria predators, Kivu predators and outgroups

i=0
for P1 in 109429 109432 103624
do
  for P2 in checkmate 130705 Ig279
  do
    for P3 in 81022 81019 81027 81063
    do
      for P4 in 64293 64394 64253
      do
        for O in 131282 Aburtoni
        do
          i=$((i+1))
          echo $i
          ./runDfoil.sh $prefix $i $P1 $P2 $P3 $P4 $O
        done
      done
    done
  done
done

# Check the introgression inference:
cut -f 32 KivInsect-KivPisc-squami-VicPisc-burtoni.*.dfoil

# Combine Dfoil results:
head -1 nubila-VicPisc-squami-KivPisc-burtoni.1.dfoil | sed 's/#//' > nubila-VicPisc-squami-KivPisc-burtoni.dfoil
for i in *.dfoil; do grep -v introgression $i >> nubila-VicPisc-squami-KivPisc-burtoni.dfoil; done
