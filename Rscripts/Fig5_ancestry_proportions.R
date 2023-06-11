require("poppr")

# Functions to get proportions of ecomorph-associated alleles found in other groups

# Read in chromosome info (to make positions additive to infer the number of genomic regions)
ref.scaff <- read.table('D:/Dropbox/victoriaGenomes/puncross.gapsEstimated.fasta.fai')
ref.scaff<-ref.scaff[,1:2]
names(ref.scaff)<-c("scaff","length")
ref.scaff$add<-c(0,cumsum(ref.scaff$length)[-length(ref.scaff$scaff)])

# Functions to extract relevant statistics

getStats<-function(ecogroup="",abb,nsp=0){
  
  if(file.exists(paste0("D:/Dropbox/victoriaGenomes/ecologicalGroups/LV.",ecogroup,".chr1-22.0.9diff.sites")))
    dat<-read.table(paste0("D:/Dropbox/victoriaGenomes/ecologicalGroups/LV.",ecogroup,".chr1-22.0.9diff.sites"),header=T)
  else{
    return("File does not exist")
  }

  # total number of SNPs with at least 0.9 allele frequency difference
  totN<-length(dat$CHROM)
  
  # Number of samples of the ecogroup
  ecoN<-max(dat[,3]/2)
  
  # Add additive positions  to infer the number of ecomorph associated regions
  dat$addPos<-dat$POS+ref.scaff[match(as.character(dat$CHROM),as.character(ref.scaff$scaff)),"add"]
  
  # spread across totChr and totRegions
  totChr<-length(unique(dat$CHROM))
  totSNPs<-nrow(dat)
  totRegions10kb<-length(which(diff(dat$addPos)>10000)+1) # regions of SNPs within 10kb
  totRegions50kb<-length(which(diff(dat$addPos)>50000)+1) # regions of SNPs within 10kb
  congo0<-length(which(dat$Frq_Congo==0&!is.na(dat$Frq_Congo)))
  propCongo0<-round(congo0/length(which(!is.na(dat$Frq_Congo))),2)
  nile0<-length(which(dat$Frq_Nile==0&!is.na(dat$Frq_Nile)))
  propNile0<-round(nile0/length(which(!is.na(dat$Frq_Nile))),2)
  parSubset<-dat[dat$N_Congo>0&dat$N_Nile>0,]
  propOnlyCongo<-round(length(which(parSubset$Frq_Congo>0&parSubset$Frq_Nile==0))/length(parSubset$CHROM),2)
  propOnlyNile<-round(length(which(parSubset$Frq_Congo==0&parSubset$Frq_Nile>0))/length(parSubset$CHROM),2)
  propBothParentals<-round(length(which(parSubset$Frq_Congo>0&parSubset$Frq_Nile>0))/length(parSubset$CHROM),2)
  propNoParentals<-round(length(which(parSubset$Frq_Congo==0&parSubset$Frq_Nile==0))/length(parSubset$CHROM),2)
  
  # Remove sites where the Southern or Northern Generalists lack data
  WestNubsDat<-dat[!is.na(dat$Frq_Victoria_insectivore_nubs)&!is.na(dat$Frq_Victoria_Northern_generalists),]
  WestNubsDat$west<-rowMeans(cbind(WestNubsDat$Frq_edward,WestNubsDat$Frq_albert,
                                   WestNubsDat$Frq_saka,WestNubsDat$Frq_kivu),na.rm=T)
  nubsOnlyProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs>0 &
                                     WestNubsDat$west==0 & WestNubsDat$Frq_Victoria_Northern_generalists==0,1])/nrow(WestNubsDat) 
  mitoOnlyProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_Northern_generalists>0 &
                                     WestNubsDat$west==0 & WestNubsDat$Frq_Victoria_insectivore_nubs==0,1])/nrow(WestNubsDat) 
  westOnlyProp<-length(WestNubsDat[WestNubsDat$west>0 & WestNubsDat$Frq_Victoria_Northern_generalists==0 &
                                     WestNubsDat$Frq_Victoria_insectivore_nubs==0,1])/nrow(WestNubsDat)
  nubsWestProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs>0 &
                                     WestNubsDat$west>0 & WestNubsDat$Frq_Victoria_Northern_generalists==0,1])/nrow(WestNubsDat)
  mitoWestProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs==0 &
                                     WestNubsDat$west>0 & WestNubsDat$Frq_Victoria_Northern_generalists>0,1])/nrow(WestNubsDat)
  mitonubsProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs>0 &
                                      WestNubsDat$west==0 & WestNubsDat$Frq_Victoria_Northern_generalists>0,1])/nrow(WestNubsDat)
  othersProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs>0 &
                                      WestNubsDat$west>0 & WestNubsDat$Frq_Victoria_Northern_generalists>0,1])/nrow(WestNubsDat)
  vicOnlyProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs==0 &
                                    WestNubsDat$west==0 & WestNubsDat$Frq_Victoria_Northern_generalists==0,1])/nrow(WestNubsDat)
  
  # Add split columns
  
  
  # Sites only found in the Congolese parental lineage
  westNubsDatCon<-WestNubsDat[WestNubsDat$Frq_Congo>0&WestNubsDat$Frq_Nile==0&(WestNubsDat$N_Congo>0&WestNubsDat$N_Nile>0),]
  
  if(nrow(westNubsDatCon)>0){
  nubsOnlyPropCon<-length(westNubsDatCon[westNubsDatCon$Frq_Victoria_insectivore_nubs>0 &
                                     westNubsDatCon$west==0 & westNubsDatCon$Frq_Victoria_Northern_generalists==0,1])/nrow(westNubsDatCon) 
  mitoOnlyPropCon<-length(westNubsDatCon[westNubsDatCon$Frq_Victoria_Northern_generalists>0 &
                                     westNubsDatCon$west==0 & westNubsDatCon$Frq_Victoria_insectivore_nubs==0,1])/nrow(westNubsDatCon) 
  westOnlyPropCon<-length(westNubsDatCon[westNubsDatCon$west>0 & westNubsDatCon$Frq_Victoria_Northern_generalists==0 &
                                     westNubsDatCon$Frq_Victoria_insectivore_nubs==0,1])/nrow(westNubsDatCon)
  nubsWestPropCon<-length(westNubsDatCon[westNubsDatCon$Frq_Victoria_insectivore_nubs>0 &
                                     westNubsDatCon$west>0 & westNubsDatCon$Frq_Victoria_Northern_generalists==0,1])/nrow(westNubsDatCon)
  mitoWestPropCon<-length(westNubsDatCon[westNubsDatCon$Frq_Victoria_insectivore_nubs==0 &
                                     westNubsDatCon$west>0 & westNubsDatCon$Frq_Victoria_Northern_generalists>0,1])/nrow(westNubsDatCon)
  mitonubsPropCon<-length(westNubsDatCon[westNubsDatCon$Frq_Victoria_insectivore_nubs>0 &
                                     westNubsDatCon$west==0 & westNubsDatCon$Frq_Victoria_Northern_generalists>0,1])/nrow(westNubsDatCon)
  othersPropCon<-length(westNubsDatCon[westNubsDatCon$Frq_Victoria_insectivore_nubs>0 &
                                   westNubsDatCon$west>0 & westNubsDatCon$Frq_Victoria_Northern_generalists>0,1])/nrow(westNubsDatCon)
  vicOnlyPropCon<-length(westNubsDatCon[westNubsDatCon$Frq_Victoria_insectivore_nubs==0 &
                                    westNubsDatCon$west==0 & westNubsDatCon$Frq_Victoria_Northern_generalists==0,1])/nrow(westNubsDatCon)
  }else
    nubsOnlyPropCon=mitoOnlyPropCon=westOnlyPropCon=nubsWestPropCon=mitoWestPropCon=mitonubsPropCon=othersPropCon=vicOnlyPropCon=0
  
  
  # Sites only found in the Upper Nile parental lineage
  westNubsDatNile<-WestNubsDat[WestNubsDat$Frq_Congo==0&WestNubsDat$Frq_Nile>0&(WestNubsDat$N_Congo>0&WestNubsDat$N_Nile>0),]
  
  if(nrow(westNubsDatNile)>0){
  nubsOnlyPropNile<-length(westNubsDatNile[westNubsDatNile$Frq_Victoria_insectivore_nubs>0 &
                                           westNubsDatNile$west==0 & westNubsDatNile$Frq_Victoria_Northern_generalists==0,1])/nrow(westNubsDatNile) 
  mitoOnlyPropNile<-length(westNubsDatNile[westNubsDatNile$Frq_Victoria_Northern_generalists>0 &
                                           westNubsDatNile$west==0 & westNubsDatNile$Frq_Victoria_insectivore_nubs==0,1])/nrow(westNubsDatNile) 
  westOnlyPropNile<-length(westNubsDatNile[westNubsDatNile$west>0 & westNubsDatNile$Frq_Victoria_Northern_generalists==0 &
                                           westNubsDatNile$Frq_Victoria_insectivore_nubs==0,1])/nrow(westNubsDatNile)
  nubsWestPropNile<-length(westNubsDatNile[westNubsDatNile$Frq_Victoria_insectivore_nubs>0 &
                                           westNubsDatNile$west>0 & westNubsDatNile$Frq_Victoria_Northern_generalists==0,1])/nrow(westNubsDatNile)
  mitoWestPropNile<-length(westNubsDatNile[westNubsDatNile$Frq_Victoria_insectivore_nubs==0 &
                                           westNubsDatNile$west>0 & westNubsDatNile$Frq_Victoria_Northern_generalists>0,1])/nrow(westNubsDatNile)
  mitonubsPropNile<-length(westNubsDatNile[westNubsDatNile$Frq_Victoria_insectivore_nubs>0 &
                                           westNubsDatNile$west==0 & westNubsDatNile$Frq_Victoria_Northern_generalists>0,1])/nrow(westNubsDatNile)
  othersPropNile<-length(westNubsDatNile[westNubsDatNile$Frq_Victoria_insectivore_nubs>0 &
                                         westNubsDatNile$west>0 & westNubsDatNile$Frq_Victoria_Northern_generalists>0,1])/nrow(westNubsDatNile)
  vicOnlyPropNile<-length(westNubsDatNile[westNubsDatNile$Frq_Victoria_insectivore_nubs==0 &
                                          westNubsDatNile$west==0 & westNubsDatNile$Frq_Victoria_Northern_generalists==0,1])/nrow(westNubsDatNile)
  }else
    nubsOnlyPropNile=mitoOnlyPropNile=westOnlyPropNile=nubsWestPropNile=mitoWestPropNile=mitonubsPropNile=othersPropNile=vicOnlyPropNile=0
  
  
  tmp<-cbind(ecogroup,totN,ecoN,totSNPs,totChr,totRegions10kb,totRegions50kb,
             congo0,propCongo0,nile0,propNile0,othersProp,
             nubsOnlyProp,westOnlyProp,nubsWestProp,vicOnlyProp,mitoOnlyProp,mitoWestProp,mitonubsProp,
             nubsOnlyPropCon,mitoOnlyPropCon,westOnlyPropCon,mitoWestPropCon,nubsWestPropCon,mitonubsPropCon,othersPropCon,vicOnlyPropCon,
             nubsOnlyPropNile,mitoOnlyPropNile,westOnlyPropNile,nubsWestPropNile,mitoWestPropNile,mitonubsPropNile,othersPropNile,vicOnlyPropNile,
             propOnlyCongo,propOnlyNile,propBothParentals,propNoParentals,abb,nsp
             )

  return(tmp)
}
getMatchFrqStats<-function(ecogroup=""){
  
  if(file.exists(paste0("D:/Dropbox/victoriaGenomes/ecologicalGroups/LV.",ecogroup,".chr1-22.matchAllFrq.sites"))){
    dat<-read.table(paste0("D:/Dropbox/victoriaGenomes/ecologicalGroups/LV.",ecogroup,".chr1-22.matchAllFrq.sites"),header=T)
  }else{
    return("File does not exist")
  }
  
  # Remove sites that are ecotype associated to get a truly independent control set
  ecoAss<-read.table(paste0("D:/Dropbox/victoriaGenomes/ecologicalGroups/LV.",ecogroup,".chr1-22.0.9diff.sites"),header=T)
  ecoSites<-paste0(ecoAss$CHROM,ecoAss$POS)
  dat<-dat[!paste0(dat$CHROM,dat$POS)%in%ecoSites,]
  
  # total number of SNPs with low allele frequency difference
  totN<-length(dat$CHROM)
  
  # Number of samples of the ecogroup
  ecoN<-max(dat[,3]/2)
  
  # Get allele frequency range
  VicMinFrq<-min(dat$Frq_Victoria,na.rm=T)
  VicMaxFrq<-max(dat$Frq_Victoria,na.rm=T)
  
  # Compute proportions of SNPs with ecomorph allele in Congo and/or Nile
  parSubset<-dat[dat$N_Congo>0&dat$N_Nile>0,]
  propOnlyCongo<-round(length(which(parSubset$Frq_Congo>0&parSubset$Frq_Nile==0))/length(parSubset$CHROM),2)
  propOnlyNile<-round(length(which(parSubset$Frq_Congo==0&parSubset$Frq_Nile>0))/length(parSubset$CHROM),2)
  propBothParentals<-round(length(which(parSubset$Frq_Congo>0&parSubset$Frq_Nile>0))/length(parSubset$CHROM),2)
  propNoParentals<-round(length(which(parSubset$Frq_Congo==0&parSubset$Frq_Nile==0))/length(parSubset$CHROM),2)
  
  
  WestNubsDat<-dat[!is.na(dat$Frq_Victoria_insectivore_nubs)&!is.na(dat$Frq_Victoria_Northern_generalists),]
  WestNubsDat$west<-rowMeans(cbind(WestNubsDat$Frq_edward,WestNubsDat$Frq_albert,
                                   WestNubsDat$Frq_saka,WestNubsDat$Frq_kivu),na.rm=T)
  nubsOnlyProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs>0 &
                                     WestNubsDat$west==0 & WestNubsDat$Frq_Victoria_Northern_generalists==0,1])/nrow(WestNubsDat) 
  mitoOnlyProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_Northern_generalists>0 &
                                     WestNubsDat$west==0 & WestNubsDat$Frq_Victoria_insectivore_nubs==0,1])/nrow(WestNubsDat) 
  westOnlyProp<-length(WestNubsDat[WestNubsDat$west>0 & WestNubsDat$Frq_Victoria_Northern_generalists==0 &
                                     WestNubsDat$Frq_Victoria_insectivore_nubs==0,1])/nrow(WestNubsDat)
  nubsWestProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs>0 &
                                     WestNubsDat$west>0 & WestNubsDat$Frq_Victoria_Northern_generalists==0,1])/nrow(WestNubsDat)
  mitoWestProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs==0 &
                                     WestNubsDat$west>0 & WestNubsDat$Frq_Victoria_Northern_generalists>0,1])/nrow(WestNubsDat)
  mitonubsProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs>0 &
                                     WestNubsDat$west==0 & WestNubsDat$Frq_Victoria_Northern_generalists>0,1])/nrow(WestNubsDat)
  othersProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs>0 &
                                   WestNubsDat$west>0 & WestNubsDat$Frq_Victoria_Northern_generalists>0,1])/nrow(WestNubsDat)
  vicOnlyProp<-length(WestNubsDat[WestNubsDat$Frq_Victoria_insectivore_nubs==0 &
                                    WestNubsDat$west==0 & WestNubsDat$Frq_Victoria_Northern_generalists==0,1])/nrow(WestNubsDat)
  

  # get significance by resampling
  westnubs<-as.numeric()
  nubs<-as.numeric()
  west<-as.numeric()
  vic<-as.numeric()
  Northern_generalists<-as.numeric()
  cong<-as.numeric()
  nile<-as.numeric()
  parentals<-as.numeric()
  lvrs<-as.numeric()
  for(i in 1:100){
    WestNubsDatSubset<-WestNubsDat[sample(1:nrow(WestNubsDat), 
                                          as.numeric(stats[stats$ecogroup==ecogroup,"totN"]), 
                                          replace=F),]
    parDatSubset<-parSubset[sample(1:length(parSubset$CHROM), 
                                   as.numeric(stats[stats$ecogroup==ecogroup,"totN"]), 
                                   replace=F),]
    westnubs<-c(westnubs,length(WestNubsDatSubset[WestNubsDatSubset$Frq_Victoria_insectivore_nubs>0 &
                                                    WestNubsDatSubset$west>0 & WestNubsDatSubset$Frq_Victoria_Northern_generalists==0,1])/nrow(WestNubsDatSubset))
    nubs<-c(nubs,length(WestNubsDatSubset[WestNubsDatSubset$Frq_Victoria_insectivore_nubs>0 &
                                                    WestNubsDatSubset$west==0 & WestNubsDatSubset$Frq_Victoria_Northern_generalists==0,1])/nrow(WestNubsDatSubset))
    west<-c(west,length(WestNubsDatSubset[WestNubsDatSubset$Frq_Victoria_insectivore_nubs==0 &
                                                    WestNubsDatSubset$west>0 & WestNubsDatSubset$Frq_Victoria_Northern_generalists==0,1])/nrow(WestNubsDatSubset))
    vic<-c(vic,length(WestNubsDatSubset[WestNubsDatSubset$Frq_Victoria_insectivore_nubs==0 &
                                            WestNubsDatSubset$west==0,1])/nrow(WestNubsDatSubset))
    Northern_generalists<-c(Northern_generalists,length(WestNubsDatSubset[WestNubsDatSubset$Frq_Victoria_insectivore_nubs==0 &
                                            WestNubsDatSubset$west==0 & WestNubsDatSubset$Frq_Victoria_Northern_generalists>0,1])/nrow(WestNubsDatSubset))
    cong<-c(cong,length(parDatSubset[parDatSubset$Frq_Congo>0 &
                                            parDatSubset$Frq_Nile==0,1])/nrow(parDatSubset))
    nile<-c(nile,length(parDatSubset[parDatSubset$Frq_Congo==0 &
                                       parDatSubset$Frq_Nile>0,1])/nrow(parDatSubset))
    parentals<-c(parentals,length(parDatSubset[parDatSubset$Frq_Congo>0 &
                                                 parDatSubset$Frq_Nile>0,1])/nrow(parDatSubset))
    lvrs<-c(lvrs,length(parDatSubset[parDatSubset$Frq_Congo==0 &
                                       parDatSubset$Frq_Nile==0,1])/nrow(parDatSubset))
  }
  westNubsSign<-stats$nubsWestProp[stats$ecogroup==ecogroup]>as.numeric(quantile(westnubs,probs = 0.95))
  westOnlySign<-stats$westOnlyProp[stats$ecogroup==ecogroup]>as.numeric(quantile(west,probs = 0.95))
  nubsOnlySign<-stats$nubsOnlyProp[stats$ecogroup==ecogroup]>as.numeric(quantile(nubs,probs = 0.95))
  vicOnlySign<-stats$vicOnlyProp[stats$ecogroup==ecogroup]>as.numeric(quantile(vic,probs = 0.95))
  Northern_generalistsOnlySign<-stats$mitoOnlyProp[stats$ecogroup==ecogroup]>as.numeric(quantile(Northern_generalists,probs = 0.95))
  
  congoSign<-stats$propOnlyCongo[stats$ecogroup==ecogroup]>as.numeric(quantile(cong,probs = 0.95))
  nileSign<-stats$propOnlyNile[stats$ecogroup==ecogroup]>as.numeric(quantile(nile,probs = 0.95))
  bothParSign<-stats$propBothParentals[stats$ecogroup==ecogroup]>as.numeric(quantile(parentals,probs = 0.95))
  noParSign<-stats$propNoParentals[stats$ecogroup==ecogroup]>as.numeric(quantile(lvrs,probs = 0.95))

  return(cbind(ecogroup,totN,ecoN,
               VicMinFrq,VicMaxFrq,nubsOnlyProp,westOnlyProp,nubsWestProp,vicOnlyProp,mitoOnlyProp,
               nubsOnlySign,westOnlySign,westNubsSign,vicOnlySign,Northern_generalistsOnlySign,
               congoSign,nileSign,bothParSign,noParSign,
               propOnlyCongo,propOnlyNile,propBothParentals,propNoParentals))
}


# Initialize dataframes
stats<-data.frame()
matchFrqStats<-data.frame()

# Enterochromis coprologus II (de)
stats<-rbind(stats,getStats("detritivore_E2","de",1))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("detritivore_E2"))

# 6 Neochromis species (ab)
stats<-rbind(stats,getStats("Neo_epilithic_algae_scraper","ab",6))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("Neo_epilithic_algae_scraper"))

# Mbipia lutea
stats<-rbind(stats,getStats("epilithic_algae_scraper_Ml","ag_M.lut",1))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("epilithic_algae_scraper_Ml"))

# Mbipia mbipi
stats<-rbind(stats,getStats("epilithic_algae_scraper_Mm","ag M.mbi",1))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("epilithic_algae_scraper_Mm"))

# Paralabidochromis short snout scraper
stats<-rbind(stats,getStats("epilithic_algae_scraper_Pshorsnsc","as",1))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("epilithic_algae_scraper_Pshorsnsc"))

# H. sp. "purple yellow"
stats<-rbind(stats,getStats("epiphytic_algae_scraper_purYel","ep",1))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("epiphytic_algae_scraper_purYel"))

# 3 Yssichromis species (plumbus, laparogramma, pyrrhocephalus)
stats<-rbind(stats,getStats("zooplanktivore_Y","zo",3))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("zooplanktivore_Y"))

# 3 Ptyochromis species (xenognathus, xenognathus rocks, fischeri)
stats<-rbind(stats,getStats("oral_sheller_Macro","os_Pty",3))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("oral_sheller_Macro"))

# P. plagiodon
stats<-rbind(stats,getStats("oral_sheller_plagiodon","os_Par",1))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("oral_sheller_plagiodon"))

# Macropleurodus bicolor
stats<-rbind(stats,getStats("oral_crusher_Macro","oc",1))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("oral_crusher_Macro"))

# 3 Labrochromis species (ismaeli, demersal MG, stone)
stats<-rbind(stats,getStats("pharyngeal.crusher","pc",3))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("pharyngeal.crusher"))

# 2 Gaurochromis species (hiatus and sp)
stats<-rbind(stats,getStats("insectivore_E1Gau","in_Gau",2))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("insectivore_E1Gau"))

# 2 Paralabidochromis species (chilotes and pseudorockpicker)
stats<-rbind(stats,getStats("insectivore_Paral","in_Par",2))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("insectivore_Paral"))

# 3 Pundamilia species (pundamilia, deepwater giant, macrocephala Makobe)
stats<-rbind(stats,getStats("insectivore_Pund","in_Pun",3))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("insectivore_Pund"))

# 2 Lithochromis (orange, yellow-chin)
stats<-rbind(stats,getStats("insectivore_Lith","in_Lit",2))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("insectivore_Lith"))

# 2 crabs eater: howesi and howesi species
stats<-rbind(stats,getStats("crabs_eater","cr",2))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("crabs_eater"))

# 5 Lipchromis species
stats<-rbind(stats,getStats("paedophage","pa",5))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("paedophage"))

# dwarf predators
stats<-rbind(stats,getStats("dwarf_predators","dw",5))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("dwarf_predators"))

# multiple piscivore species
stats<-rbind(stats,getStats("piscivore","pi",16))
matchFrqStats<-rbind(matchFrqStats,getMatchFrqStats("piscivore"))

# Combine stats and significance info for the same ecomorphs
stats<-cbind(stats,matchFrqStats[,grepl(names(matchFrqStats),pattern="Sign")])

# Add abbreviations as rownames
rownames(stats)<-stats$abb


# Proportion of sites that the allele may have evolved de novo in Lake Victoria
stats$vicOnlyProp


# png(file="D:/Dropbox/victoriaGenomes/ecologicalGroups/ecomorph_categories_simplified.png",
#     width=1300,height=700)
pdf(file="D:/Dropbox/victoriaGenomes/ecologicalGroups/ecomorph_categories_simplified.pdf",
    width=13,height=7)
par(mfrow=c(3,1),mar=c(2,4,1,1),xaxs = "i",mgp=c(1.9,0.5,0))
bp<-barplot(matrix(rbind(as.double(stats[,"nubsOnlyProp"]),
                         as.double(stats[,"mitoOnlyProp"]),
                         as.double(stats[,"westOnlyProp"]),
                         as.double(stats[,"mitonubsProp"]),
                         as.double(stats[,"nubsWestProp"]),
                         as.double(stats[,"mitoWestProp"]),
                     as.double(stats[,"othersProp"]),
                     as.double(stats[,"vicOnlyProp"])),
            nrow = 8),names.arg = stats$abb,
        ylim=c(0,1.6), col=c("orange3","orange","cornflowerblue","grey30","grey50","grey70","grey90","white"),
        ylab="Proportion of SNPs                   ",yaxt="n")
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1),las=2)
legend("top",ncol=3,fill=c("orange3","orange","cornflowerblue","grey30","grey50","grey70","grey90","white"),
       legend=c("A. nubila only","Northern generalists only","Western Lakes only","A. nubila and northern generalists",
                "nubila and Western","N generalists and Western","nubila,N generalists,Western","only in Victoria"))
text(x=rep(bp,each=3),y=as.double(apply(matrix(rbind(as.double(stats[,"nubsOnlyProp"]),
                                                     as.double(stats[,"mitoOnlyProp"]),
                                                     as.double(stats[,"westOnlyProp"])),
                                               nrow = 3),2,cumsum))-0.02,
     labels = ifelse(as.logical(as.matrix(t(cbind(
       stats$nubsOnlySign,stats$Northern_generalistsOnlySign,stats$westOnlySign))))==TRUE,"*",""),cex=2
)
text(bp,y=1.05,cex=0.8,labels = paste0(stats$ecoN,"/",stats$nsp))
text(0.1,y=1.4,cex=0.8,labels = "individuals N/species N",xpd=T)
text(bp,y=1.15,cex=0.8,labels = paste0(stats$totSNPs,"/",stats$totRegions10kb))
text(0.1,y=1.5,cex=0.8,labels = "SNPs/regions",xpd=T)
arrows(x0 = 0,x1=bp[1],y0=1.4,y1=1.2,length = 0.15)

# Upper Nile / Congo allele presence
bp<-barplot(matrix(rbind(as.double(stats[,"propOnlyCongo"]),
                         as.double(stats[,"propOnlyNile"]),
                         as.double(stats[,"propBothParentals"]),
                         as.double(stats[,"propNoParentals"])),
                   nrow = 4),names.arg = stats$abb,
            ylim=c(0,1.6), col=c("red","blue","grey30","white"),
            ylab="Proportion of SNPs                  ",yaxt="n")
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1),las=2)
legend("top",fill=c("red","blue","grey30","white"),
       legend=c("Congolese lineage only","Upper Nile lineage only",
                "both Congo-Nilotic parentals","absent in Congo-Nilotic parentals"),ncol=2)
text(x=rep(bp,each=2),y=as.double(apply(matrix(rbind(as.double(stats[,"propOnlyCongo"]),
                                                     as.double(stats[,"propOnlyNile"])),
                                               nrow = 2),2,cumsum))-0.02,
     labels = ifelse(as.logical(as.matrix(t(cbind(
       matchFrqStats$congoSign ,matchFrqStats$nileSign))))==TRUE,"*",""),cex=2)

# Split barplots
bp<-barplot(matrix(rbind(as.double(stats[,"nubsOnlyPropCon"]),
                         as.double(stats[,"mitoOnlyPropCon"]),
                         as.double(stats[,"westOnlyPropCon"]),
                         as.double(stats[,"mitonubsPropCon"]),
                         as.double(stats[,"nubsWestPropCon"]),
                         as.double(stats[,"mitoWestPropCon"]),
                         as.double(stats[,"othersPropCon"]),
                         as.double(stats[,"vicOnlyPropCon"]),
                         as.double(stats[,"nubsOnlyPropNile"]),
                         as.double(stats[,"mitoOnlyPropNile"]),
                         as.double(stats[,"westOnlyPropNile"]),
                         as.double(stats[,"mitonubsPropNile"]),
                         as.double(stats[,"nubsWestPropNile"]),
                         as.double(stats[,"mitoWestPropNile"]),
                         as.double(stats[,"othersPropNile"]),
                         as.double(stats[,"vicOnlyPropNile"])),
                   nrow = 8),names.arg = "",
            ylim=c(0,1), col=c("orange3","orange","cornflowerblue","grey30","grey50","grey70","grey90","white"),
            ylab="Proportion of SNPs                   ",yaxt="n", space=c(0.2,0))
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels = c(0,0.2,0.4,0.6,0.8,1),las=2)
axis(1,at=c(1.2,3.4,5.6,7.8,10.0,12.2,14.4,16.6,18.8,21.0,23.2,25.4,27.6,29.8,32.0,34.2,36.4,38.6,40.8),labels=stats$abb)


dev.off()

# Save the data as spread sheet
write.csv(stats,file="D:/Dropbox/victoriaGenomes/ecologicalGroups/Fig5_ancestry_proportions.csv")
