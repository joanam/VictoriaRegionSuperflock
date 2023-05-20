
### Fig. S10: Piscivore-associated sites #####

frq<-read.table("D:/Dropbox/victoriaGenomes/ecologicalGroups/LV.piscivore.chr1-22.0.9diff.sites",header=T)


# difference between piscivores and others in Albert/Kivu
frq$diffAlb<-frq$Frq_Albert_piscivore-frq$Frq_Albert_nonpiscivore
frq$diffKiv<-frq$Frq_Kivu_piscivore-frq$Frq_Kivu_nonpiscivore

# Get matching allele frequencies 
ecogroup="piscivore"
matchDat<-read.table(paste0("D:/Dropbox/victoriaGenomes/ecologicalGroups/LV.",ecogroup,".chr1-22.matchAllFrq.sites"),header=T)

# Remove sites that are ecotype associated to get a truly independent control set
ecoSites<-paste0(frq$CHROM,frq$POS)
matchDat<-matchDat[!paste0(matchDat$CHROM,matchDat$POS)%in%ecoSites,]

# Get upper threshold for expected difference between piscivores and non-piscivores in Edward/Kivu
matchDat$diffAlb<-matchDat$Frq_Albert_piscivore-matchDat$Frq_Albert_nonpiscivore
matchDat$diffKiv<-matchDat$Frq_Kivu_piscivore-matchDat$Frq_Kivu_nonpiscivore
threshAlb<-quantile(matchDat$diffAlb[matchDat$N_Albert_piscivore>0],probs = 0.95,na.rm=T)
threshKiv<-quantile(matchDat$diffKiv[matchDat$N_Kivu_piscivore>2],probs = 0.95,na.rm=T)



### Get genes ####
{
  #Get all genes
  mRNA<-read.csv("D:/Dropbox/refGenomes/NCBI_tilapia_nyererei_mRNA.csv",header=T,sep=";")
  mRNA<-mRNA[!grepl(mRNA$CHR,pattern="scaffold"),]
  
  # Get genes coinciding with piscivore-associated regions
  genesPisc<-data.frame()
  for(i in 1:length(frq$CHROM)){
    win<-frq[i,1:2]
    start<-as.integer(as.character(win$POS))-25000
    end<-as.integer(as.character(win$POS))+25000
    chr<-as.character(win$CHROM)
    chrGenes<-mRNA[mRNA$CHR==chr,]
    chrGenes<-chrGenes[(chrGenes$start<=start & chrGenes$end>=start) |
                         (chrGenes$start<=end & chrGenes$end>=end) |
                         (chrGenes$start>=start & chrGenes$end<=end) |
                         (chrGenes$start<=start & chrGenes$end>=end),]
    
    # If more than one genes occur in the interval, get the closest one
    if(length(chrGenes$CHR)>1){
      dist<-apply(MARGIN = 1,X = cbind(abs(chrGenes$start-win$POS),abs(chrGenes$end-win$POS)),FUN = min)
      chrGenes<-chrGenes[which(dist==min(dist)),]
    }
    # Add the gene to the gene collection
    if(length(chrGenes$CHR)>0) genesPisc<-rbind(genesPisc,cbind(chrGenes,win))
  }
  genesSNPs<-aggregate(x = genesPisc['POS'], by=genesPisc[c("CHROM","start","end")], toString)
  names(genesSNPs)[4]<-"SNP_POS"
  genesSNPs$SNP_POS
  genesPisc<-unique(genesPisc[,1:15])
  genesPisc<-merge(genesSNPs,genesPisc,by=c("CHROM","start","end"))
  genesPisc<-as.data.frame(genesPisc)
  
  require(gtools)
  genesPisc<-genesPisc[order(genesPisc$start),]
  genesPisc<-genesPisc[mixedorder(genesPisc$CHR),]
  
  genesPisc<-read.table("D:/Dropbox/victoriaGenomes/ecologicalGroups/piscivore-genes.txt",header=T,sep="\t")
}


# Add additive positions only for chromosomes that have sites
ref.scaff <- read.table('D:/Dropbox/CichlidGenomesProject/refGenome/puncross.gapsEstimated.fasta.fai')
ref.scaff<-ref.scaff[,1:2]
names(ref.scaff)<-c("scaff","length")
ref.scaff<-ref.scaff[grepl(ref.scaff$scaff,pattern="chr"),]
ref.scaff<-ref.scaff[ref.scaff$scaff!="chrM",]
refSubset<-ref.scaff[ref.scaff$scaff%in%frq$CHROM,]
refSubset$add<-c(0,cumsum(refSubset$length)[-length(refSubset$scaff)])

# Get additive positions
frq$pos_add<-frq$POS+refSubset[match(as.character(frq$CHROM),as.character(refSubset$scaff)),"add"]

# Plot allele frequencies for Fig. S10
png(file="D:/Dropbox/victoriaGenomes/ecologicalGroups/piscivore-associated-SNPs.png",
     width=31,height = 19, units = 'cm', res = 300)
par(mfrow=c(5,1),oma=c(5,0,4,0),mar=c(0,5,1,0),xaxs="i")

{

# Lake Victoria
plot(frq$pos_add,frq$Frq_Victoria_piscivore,ylim=c(0,1),las=2,
     xaxt="n",xlab="",ylab="Allele frequency",
     xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]),col="darkorange",pch=19,cex=1.5)
points(frq$pos_add,frq$Frq_Victoria_non_piscivore,pch=8)
abline(v=refSubset$add)
legend(x=52000000,y=1,legend=c("25 piscivores","248 nonpiscivores"),
       pch=c(19,8),box.lty=0,bg="white",col=c("darkorange","black"),title="  Victoria Radiation", title.adj=0)
text(labels = "A",x=-9500000,y=1,xpd=T,cex=1.3)

# Get additive positions of relevant genes
genesStart<-genesPisc$start+
  refSubset[match(as.character(genesPisc$CHR),as.character(refSubset$scaff)),"add"]
genesEnd<-genesPisc$end+
  refSubset[match(as.character(genesPisc$CHR),as.character(refSubset$scaff)),"add"]

# Add gene names to the plot
text(seq(grconvertX(0.05,"npc","user"),
         grconvertX(0.9,"npc","user"),
         length.out=dim(genesPisc)[1]),
     grconvertY(1.03,"nic","user"),
     genesPisc$short_name,cex=1.3,xpd=NA,
     adj=0.5,srt=0)

# Link gene names to genes
segments(x0 = genesStart+(genesEnd-genesStart)/2,
         y0=grconvertY(0.98,"nic","user"),
         x1=seq(grconvertX(0.05,"npc","user"),
                grconvertX(0.9,"npc","user"),
                length.out=dim(genesPisc)[1]),
         y1=grconvertY(1.02,"nic","user"),
         col="black",
         xpd=NA,lwd=1)

# A nubila
plot(frq$pos_add,frq$Frq_Victoria_insectivore_nubs,ylim=c(0,1),pch=8,
     xaxt="n",xlab="",ylab="Allele frequency",las=2,
     xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]))
abline(v=refSubset$add)
legend(x=52000000,y=1,pch=8,legend=c("11 A. nubila                        "),
       col=c("black"),box.lty=0,bg="white",title="  streams South of Lake Victoria", title.adj=0)
text(labels = "B",x=-9500000,y=1,xpd=T,cex=1.3)

# Northern Generalists
plot(frq$pos_add,frq$Frq_Victoria_mitoXysti,ylim=c(0,1),pch=8,
     xaxt="n",xlab="",ylab="Allele frequency",las=2,
     xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]))
legend(x=52000000,y=1,pch=8,legend=c("4 individuals       "),
       box.lty=0,bg="white",title="  Northern Generalists", title.adj=0)
abline(v=refSubset$add)
text(labels = "C",x=-9500000,y=1,xpd=T,cex=1.3)

# Kivu
plot(frq$pos_add,frq$diffKiv,
     col=ifelse(frq$N_Kivu_piscivore<3,NA,ifelse(frq$diffKiv>=threshKiv,"black","grey")),xaxt="n",xlab="",ylab="Freq. difference",
     las=2,pch=19,xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]))
legend(x=52000000,y=1,legend=c("frequency difference >= 0.47","frequency difference < 0.47                     "),
       pch=c(19),col=c("black","grey"),box.lty=0,bg="white",
       title="  Lake Kivu (3 piscivores, 41 others)", title.adj=0)
abline(v=refSubset$add)
abline(h=0)
text(labels = "D",x=-9500000,y=1,xpd=T,cex=1.3)

# Albert
plot(frq$pos_add,frq$diffAlb,
     col=ifelse(frq$N_Albert_piscivore<2,NA,ifelse(frq$diffAlb>=threshAlb,"black","grey")),xaxt="n",xlab="",ylab="Freq. difference",
     las=2,pch=19,xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]),ylim=c(-1,1))
legend(x=52000000,y=1,legend=c("frequency difference >= 0.5","frequency difference < 0.5                     "),
       pch=c(19),col=c("black","grey"),box.lty=0,bg="white",
       title="  Lake Albert (1 piscivore, 10 others)", title.adj=0)
abline(h=0)
abline(v=refSubset$add)
text(labels = "E",x=-9500000,y=1,xpd=T,cex=1.3)
axis(1,frq$pos_add,labels=round(frq$POS/1000000,digits=2))
axis(1,at=refSubset$add+refSubset$length/2,refSubset$scaff,tick=F,line=1.3,cex.axis=1.3)
mtext(1,text = "Chromosome position [Mb]",line=4)
}

dev.off()


### Fig. S11: Paedophages ####

frq<-read.table("D:/Dropbox/victoriaGenomes/ecologicalGroups/LV.paedophage.chr1-22.0.9diff.sites",header=T)

# difference between paedophage and others in Edward/Kivu
frq$diff<-frq$Frq_EdKi_paedophage-frq$Frq_EdKi_nonpaedophage

# Add additive positions  to infer the number of ecomorph associated regions
# Read in chromosome info (to make positions additive to infer the number of genomic regions)
ref.scaff <- read.table('D:/Dropbox/CichlidGenomesProject/refGenome/puncross.gapsEstimated.fasta.fai')
ref.scaff<-ref.scaff[,1:2]
names(ref.scaff)<-c("scaff","length")
ref.scaff$add<-c(0,cumsum(ref.scaff$length)[-length(ref.scaff$scaff)])

frq$addPos<-frq$POS+ref.scaff[match(as.character(frq$CHROM),as.character(ref.scaff$scaff)),"add"]

# spread across totChr and totRegions
totChr<-length(unique(frq$CHROM))
totSNPs<-nrow(frq)
totRegions10kb<-length(which(diff(frq$addPos)>10000)+1) # regions of SNPs within 10kb


# Get matching allele frequencies 
ecogroup="paedophage"
matchDat<-read.table(paste0("D:/Dropbox/victoriaGenomes/ecologicalGroups/LV.",ecogroup,".chr1-22.matchAllFrq.sites"),header=T)

# Remove sites that are ecotype associated to get a truly independent control set
ecoSites<-paste0(frq$CHROM,frq$POS)
matchDat<-matchDat[!paste0(matchDat$CHROM,matchDat$POS)%in%ecoSites,]

# Get upper threshold for expected difference between paedophages and non-paedophages in Edward/Kivu
matchDat<-matchDat[matchDat$N_EdKi_paedophage>3,]
matchDat$diff<-matchDat$Frq_EdKi_paedophage-matchDat$Frq_EdKi_nonpaedophage
thresh<-quantile(matchDat$diff,probs = 0.95,na.rm=T)

# Get genes coinciding with paedophage-associated regions
genesPaedo<-data.frame()
frqDiffEdKi<-frq[frq$diff>0.3,]
for(i in 1:length(frqDiffEdKi$CHROM)){
  win<-frqDiffEdKi[i,1:2]
  start<-as.integer(as.character(win$POS))-25000
  end<-as.integer(as.character(win$POS))+25000
  chr<-as.character(win$CHROM)
  chrGenes<-mRNA[mRNA$CHR==chr,]
  chrGenes<-chrGenes[(chrGenes$start<=start & chrGenes$end>=start) |
                       (chrGenes$start<=end & chrGenes$end>=end) |
                       (chrGenes$start>=start & chrGenes$end<=end) |
                       (chrGenes$start<=start & chrGenes$end>=end),]
  
  # If more than one genes occur in the interval, get the closest one
  if(length(chrGenes$CHR)>1){
    dist<-apply(MARGIN = 1,X = cbind(abs(chrGenes$start-win$POS),abs(chrGenes$end-win$POS)),FUN = min)
    chrGenes<-chrGenes[which(dist==min(dist)),]
  }
  # Add the gene to the gene collection
  if(length(chrGenes$CHR)>0) genesPaedo<-rbind(genesPaedo,cbind(chrGenes,win))
}
genesSNPs<-aggregate(x = genesPaedo['POS'], by=genesPaedo[c("CHROM","start","end")], toString)
names(genesSNPs)[4]<-"SNP_POS"
genesSNPs$SNP_POS
genesPaedo<-unique(genesPaedo[,1:15])
genesPaedo<-merge(genesSNPs,genesPaedo,by=c("CHROM","start","end"))
genesPaedo<-as.data.frame(genesPaedo)

require(gtools)
genesPaedo<-genesPaedo[order(genesPaedo$start),]
genesPaedo<-genesPaedo[mixedorder(genesPaedo$CHR),]


# png(file="D:/Dropbox/victoriaGenomes/ecologicalGroups/paedophage-associated-SNPs.png",
#     width=31,height = 17, units = 'cm', res = 300)
par(mar=c(0,4,1,0),mfrow=c(4,1),oma=c(6,3,10,1),xaxs="i")

{
  # Add additive positions only for chromosomes that have sites
  ref.scaff <- read.table('D:/Dropbox/CichlidGenomesProject/refGenome/puncross.gapsEstimated.fasta.fai')
  ref.scaff<-ref.scaff[,1:2]
  names(ref.scaff)<-c("scaff","length")
  ref.scaff<-ref.scaff[grepl(ref.scaff$scaff,pattern="chr"),]
  ref.scaff<-ref.scaff[ref.scaff$scaff!="chrM",]
  refSubset<-ref.scaff[ref.scaff$scaff%in%frq$CHROM,]
  refSubset$add<-c(0,cumsum(refSubset$length)[-length(refSubset$scaff)])
  
  # Get additive positions
  frq$pos_add<-frq$POS+refSubset[match(as.character(frq$CHROM),as.character(refSubset$scaff)),"add"]
  
  # Lake Victoria
  plot(frq$pos_add,frq$Frq_Victoria_paedophage,ylim=c(0,1),las=2,
       xaxt="n",xlab="",ylab="Allele frequency",
       xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]),col="darkorange",pch=19,cex=1.5)
  points(frq$pos_add,frq$Frq_Victoria_non_paedophage,pch=8)
  abline(v=refSubset$add)
  legend(0,0.9,legend=c("10 paedophages","338 nonpaedophages"),
         pch=c(19,8),col=c("darkorange","black"),box.lty=0,bg="white",title="  Victoria", title.adj=0)
  mtext(2,text = "A",las=2,line=3,at = 1,cex=1.1)
  
  
  # # Get additive positions of relevant genes
  genesStart<-genesPaedo$start+
    refSubset[match(as.character(genesPaedo$CHR),as.character(refSubset$scaff)),"add"]
  genesEnd<-genesPaedo$end+
    refSubset[match(as.character(genesPaedo$CHR),as.character(refSubset$scaff)),"add"]

  # Add gene names to the plot
  text(x = seq(grconvertX(0,"npc","user"),
           grconvertX(0.9,"npc","user"),
           length.out=dim(genesPaedo)[1]),
      y=grconvertY(1.01,"nic","user"),
       genesPaedo$gene,cex=1,xpd=NA,
       adj=0,srt=20)

  # Link gene names to genes
  segments(x0 = genesStart+(genesEnd-genesStart)/2,
           y0=grconvertY(0.98,"nic","user"),
           x1=seq(grconvertX(0,"npc","user"),
                  grconvertX(0.9,"npc","user"),
                  length.out=dim(genesPaedo)[1]),
           y1=grconvertY(1.01,"nic","user"),
           col="black",
           xpd=NA,lwd=0.2)
  
  # A nubila
  plot(frq$pos_add,frq$Frq_Victoria_insectivore_nubs,ylim=c(0,1),pch=8,
       xaxt="n",xlab="",ylab="Allele frequency",las=2,
       xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]))
  abline(v=refSubset$add)
  legend("topleft",pch=8,legend=c("11 Astatotilapia nubila                                      "),
         box.lty=0,bg="white",title="  swamps/streams South of Lake Victoria", title.adj=0)
  mtext(2,text = "B",las=2,line=3,at = 1,cex=1.1)
  
  # Northern Generalists
  plot(frq$pos_add,frq$Frq_Victoria_mitoXysti,ylim=c(0,1),pch=8,
       xaxt="n",xlab="",ylab="Allele frequency",las=2,
       xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]))
  abline(v=refSubset$add)
  legend(x="topleft",pch=8,legend=c("4 individuals                "),
         box.lty=0,bg="white",title="  Northern Generalists", 
         title.adj=0)
  mtext(2,text = "C",las=2,line=3,at = 1,cex=1.1)
  
  
  # Kivu+Edward
  plot(frq$pos_add,frq$diff,
       col=ifelse(frq$N_EdKi_paedophage<3,NA,ifelse(frq$diff>thresh,"black","grey")),xaxt="n",xlab="",ylab="Frequency difference",
       las=2,pch=19,xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]))
  abline(v=refSubset$add)
  abline(h=0)
  legend("topleft",legend=c("frequency difference > 0.28","frequency difference <= 0.28                     "),
         pch=c(19),col=c("black","grey"),box.lty=0,bg="white",
         title="  Lake Kivu+Edward (2 paedophages, 60 others)", title.adj=0)
  mtext(2,text = "D",las=2,line=3,at = 1,cex=1.1)
  
  axis(1,at=refSubset$add+refSubset$length/2,refSubset$scaff,tick=F,line=1.3,cex.axis=1.3)
  axis(1,frq$pos_add,labels=round(frq$POS/1000000,digits=2))
  mtext(1,text = "Chromosome position [Mb]",line=4)
  
}
  
dev.off()


### Fig. S12: Neochromis ####
frq<-read.table("D:/Dropbox/victoriaGenomes/ecologicalGroups/LV.Neo_epilithic_algae_scraper.chr1-22.0.9diff.sites",header=T)

# difference between Neochromis and others in Albert/Kivu
frq$diffAlb<-frq$Frq_Albert_epilithic_algae_scraper-frq$Frq_Albert_non_epilithic_algae_scraper
frq$diffKiv<-frq$Frq_Kivu_epilithic_algae_scraper-frq$Frq_Kivu_non_epilithic_algae_scraper

# Get matching allele frequencies 
ecogroup="epilithic_algae_scraper"
matchDat<-read.table(paste0("D:/Dropbox/victoriaGenomes/ecologicalGroups/LV.",ecogroup,".chr1-22.matchAllFrq.sites"),header=T)

# Get upper threshold for expected difference between paedophages and non-paedophages in Edward/Kivu
matchDat$diffAlb<-matchDat$Frq_Albert_epilithic_algae_scraper-matchDat$Frq_Albert_non_epilithic_algae_scraper
matchDat$diffKiv<-matchDat$Frq_Kivu_epilithic_algae_scraper-matchDat$Frq_Kivu_non_epilithic_algae_scraper
threshAlb<-quantile(matchDat$diffAlb[matchDat$N_Albert_epilithic_algae_scraper>0],probs = 0.95,na.rm=T)
threshKiv<-quantile(matchDat$diffKiv[matchDat$N_Kivu_epilithic_algae_scraper>2],probs = 0.95,na.rm=T)


#png(file="D:/Dropbox/victoriaGenomes/ecologicalGroups/Neochromis-associated-SNPs.png",
#    width=31,height = 17, units = 'cm', res = 300)
layout(matrix(1:10,nrow=5))
par(mar=c(0,6,1,0),oma=c(5,4,10,1),xaxs="i")

{
  # Add additive positions only for chromosomes that have sites
  ref.scaff <- read.table('D:/Dropbox/CichlidGenomesProject/refGenome/puncross.gapsEstimated.fasta.fai')
  ref.scaff<-ref.scaff[,1:2]
  names(ref.scaff)<-c("scaff","length")
  ref.scaff<-ref.scaff[grepl(ref.scaff$scaff,pattern="chr"),]
  ref.scaff<-ref.scaff[ref.scaff$scaff!="chrM",]
  refSubset<-ref.scaff[ref.scaff$scaff%in%frq$CHROM,]
  refSubset$add<-c(0,cumsum(refSubset$length)[-nrow(refSubset)])
  refSubset$add<-refSubset$add+c(0,cumsum(rep(2000000,times=nrow(refSubset)-1)))
  
  # Get additive positions
  frq$pos_add<-frq$POS+refSubset[match(as.character(frq$CHROM),as.character(refSubset$scaff)),"add"]
  
  # Get genes coinciding with Neochromis-associated regions
  genesNeo<-data.frame()
  for(i in 1:length(frq$CHROM)){
    win<-frq[i,1:2]
    start<-as.integer(as.character(win$POS))-25000
    end<-as.integer(as.character(win$POS))+25000
    chr<-as.character(win$CHROM)
    chrGenes<-mRNA[mRNA$CHR==chr,]
    chrGenes<-chrGenes[(chrGenes$start<=start & chrGenes$end>=start) |
                         (chrGenes$start<=end & chrGenes$end>=end) |
                         (chrGenes$start>=start & chrGenes$end<=end) |
                         (chrGenes$start<=start & chrGenes$end>=end),]
    
    # If more than one gene occurs in the interval, get the closest one
    if(length(chrGenes$CHR)>1){
      dist<-apply(MARGIN = 1,X = cbind(abs(chrGenes$start-win$POS),abs(chrGenes$end-win$POS)),FUN = min)
      chrGenes<-chrGenes[which(dist==min(dist)),]
    }
    # Add the gene to the gene collection
    if(length(chrGenes$CHR)>0) genesNeo<-rbind(genesNeo,cbind(chrGenes,win))
  }
  genesSNPs<-aggregate(x = genesNeo['POS'], by=genesNeo[c("CHROM","start","end")], toString)
  names(genesSNPs)[4]<-"SNP_POS"
  genesSNPs$SNP_POS
  genesNeo<-unique(genesNeo[,1:15])
  genesNeo<-merge(genesSNPs,genesNeo,by=c("CHROM","start","end"))
  genesNeo<-as.data.frame(genesNeo)
  
  require(gtools)
  genesNeo<-genesNeo[order(genesNeo$start),]
  genesNeo<-genesNeo[mixedorder(genesNeo$CHR),]
  
  
  # Lake Victoria
  plot(frq$pos_add,frq$Frq_Victoria_epilithic_algae_scraper_Neo,ylim=c(0,1),las=2,
       xaxt="n",xlab="",ylab="Allele frequency",
       xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]),col="darkorange",pch=19,cex=1.5)
  points(frq$pos_add,frq$Frq_Victoria_non_epilithic_algae_scraper,pch=8)
  abline(v=c(refSubset$add,refSubset$add+refSubset$length))
  legend("topleft",legend=c("12 algae browser","215 others"),
         pch=c(19,8),col=c("darkorange","black"),bty="n",title="  Victoria", title.adj=0)
  mtext(2,text = "A",las=2,line=3,at = 1,cex=1.1)
  
  
  # Get additive positions of relevant genes
  genesStart<-genesNeo$start+
    refSubset[match(as.character(genesNeo$CHR),as.character(refSubset$scaff)),"add"]
  genesEnd<-genesNeo$end+
    refSubset[match(as.character(genesNeo$CHR),as.character(refSubset$scaff)),"add"]
  
  # Add gene names to the plot
  text(x = seq(grconvertX(0.2,"npc","user"),
               grconvertX(0.8,"npc","user"),
               length.out=dim(genesNeo)[1]),
       y=grconvertY(1.01,"nic","user"),
       genesNeo$gene,cex=0.8,xpd=NA,
       adj=0,srt=20)
  
  # Link gene names to genes
  segments(x0 = genesStart+(genesEnd-genesStart)/2,
           y0=grconvertY(0.975,"nic","user"),
           x1=seq(grconvertX(0.2,"npc","user"),
                  grconvertX(0.8,"npc","user"),
                  length.out=dim(genesNeo)[1]),
           y1=grconvertY(1.0,"nic","user"),
           xpd=NA,lwd=0.2)
  
  # A nubila
  plot(frq$pos_add,frq$Frq_Victoria_insectivore_nubs,ylim=c(0,1),pch=8,
       xaxt="n",xlab="",ylab="Allele frequency",las=2,
       xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]))
  abline(v=c(refSubset$add,refSubset$add+refSubset$length))
  legend("topleft",pch=8,legend=c("11 Astatotilapia nubila                                      "),
         bty="n",title="  swamps/streams South of Lake Victoria", title.adj=0)
  mtext(2,text = "B",las=2,line=3,at = 1,cex=1.1)
 
  # Northern Generalists
  plot(frq$pos_add,frq$Frq_Victoria_Northern_generalists,ylim=c(0,1),pch=8,
       xaxt="n",xlab="",ylab="Allele frequency",las=2,
       xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]))
  abline(v=c(refSubset$add,refSubset$add+refSubset$length))
  legend("topleft",pch=8,legend=c("4 Northern Generalists  "),
         bty="n",title="  Northern Generalists", title.adj=0)
  mtext(2,text = "C",las=2,line=3,at = 1,cex=1.1)
  
  # Albert
  plot(frq$pos_add,frq$diffAlb,ylim=c(-0.1,1),pch=19,cex=1.5,
       xaxt="n",xlab="",las=2,
       col=ifelse(frq$N_Albert_epilithic_algae_scraper<2,NA,ifelse(frq$diffAlb>=threshAlb,"black","grey")),ylab="Frequency difference",
       xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]))
  abline(v=c(refSubset$add,refSubset$add+refSubset$length))
  abline(h=0)
  legend(x=1,y=1,bty="n",legend=c("frequency difference >= 0.44","frequency difference < 0.44                     "),
         pch=c(19),col=c("black","grey"),box.lty=0,bg="white",
         title="  Lake Albert (1 epilithic algae scraper, 10 others)", title.adj=0)
  mtext(2,text = "D",las=2,line=3,at = 1,cex=1.1)
  

  # Kivu
  plot(frq$pos_add,frq$diffKiv,ylim=c(-0.1,1),pch=19,cex=1.5,
       xaxt="n",xlab="",las=2,
       col=ifelse(frq$N_Kivu_epilithic_algae_scraper<6,NA,ifelse(frq$diffKiv>=threshKiv,"black","grey")),ylab="Frequency difference",
       xlim=c(0,refSubset$add[length(refSubset$add)]+refSubset$length[length(refSubset$add)]))
  abline(v=c(refSubset$add,refSubset$add+refSubset$length))
  abline(h=0)
  legend(x=1,y=1,legend=c("frequency difference >= 0.19","frequency difference < 0.19                     "),
         pch=c(19),bty="n",col=c("black","grey"),box.lty=0,bg="white",
         title="  Lake Kivu (4 epilithic algae scrapers, 40 others)", title.adj=0)
  mtext(2,text = "E",las=2,line=3,at = 1,cex=1.1)
  axis(1,at=refSubset$add+refSubset$length/2,refSubset$scaff,tick=F,line=1.3)
  axis(1,frq$pos_add,labels=round(frq$POS/1000000,digits=2))
  mtext(1,text = "Chromosome position [Mb]",line=4)
}

# Zoom into chr11
{
  frqZoom<-frq0.8neo[frq0.8neo$CHROM=="chr11",]
  genesNeo<-genesNeo[genesNeo$CHROM=="chr11",]

  
  # Lake Victoria
  plot(frqZoom$POS,frqZoom$Frq_Victoria_epilithic_algae_scraper_Neo,ylim=c(0,1),las=2,
       xaxt="n",xlab="",ylab="Allele frequency",xlim=c(0,111000),col="darkorange",pch=19,cex=1.5)
  points(frqZoom$POS,frqZoom$Frq_Victoria_non_epilithic_algae_scraper,pch=8)
  abline(v=refSubset$add)
  legend("topleft",legend=c("12 Neochromis","215 non-epilithic algae scrapers"),
         pch=c(19,8),col=c("darkorange","black"),bty="n",title="  Victoria", title.adj=0)
  mtext(2,text = "E",las=2,line=3,at = 1,cex=1.1)
  
  
  # Get additive positions of relevant genes
  genesStart<-genesNeo$start
  genesEnd<-genesNeo$end
  
  # Add gene names to the plot
  text(x = seq(grconvertX(0,"npc","user"),
               grconvertX(0.8,"npc","user"),
               length.out=dim(genesNeo)[1]),
       y=grconvertY(1.01,"nic","user"),
       genesNeo$gene,cex=0.8,xpd=NA,
       adj=0,srt=0)
  
  # Link gene names to genes
  segments(x0 = genesStart+(genesEnd-genesStart)/2,
           y0=grconvertY(0.975,"nic","user"),
           x1=seq(grconvertX(0.2,"npc","user"),
                  grconvertX(0.8,"npc","user"),
                  length.out=dim(genesNeo)[1]),
           y1=grconvertY(1.0,"nic","user"),
           xpd=NA,lwd=0.2)
  
  # A nubila
  plot(frqZoom$POS,frqZoom$Frq_Victoria_insectivore_nubs,ylim=c(0,1),pch=8,
       xaxt="n",xlab="",ylab="Allele frequency",las=2,xlim=c(0,111000))
  abline(v=refSubset$add)
  legend("topleft",pch=8,legend=c("11 Astatotilapia nubila                                      "),
         bty="n",title="  swamps/streams South of Lake Victoria", title.adj=0)
  mtext(2,text = "F",las=2,line=3,at = 1,cex=1.1)
  
  # Northern Generalists
  plot(frqZoom$POS,frqZoom$Frq_Victoria_Northern_generalists,ylim=c(0,1),pch=8,
       xaxt="n",xlab="",ylab="Allele frequency",las=2,xlim=c(0,111000))
  abline(v=refSubset$add)
  legend("topleft",pch=8,legend=c("4 Northern Generalists"),
         bty="n",title="  Northern Generalists", title.adj=0)
  mtext(2,text = "G",las=2,line=3,at = 1,cex=1.1)
  
  # # Saka
  # plot(frqZoom$POS,frqZoom$Frq_saka,ylim=c(0,1),pch=19,cex=1.5,
  #      xaxt="n",xlab="",las=2,xlim=c(0,110000))
  # abline(v=refSubset$add)
  # abline(h=0)
  # legend(x=1,y=1,legend=c("11 epilithic algae scraper"),
  #        pch=c(19),col=c("black"),box.lty=0,bg="white",
  #        title="  Lake Saka", title.adj=0)
  # mtext(2,text = "D",las=2,line=3,at = 1,cex=1.1)
  # 
  
  # Albert
  plot(frqZoom$POS,frqZoom$Frq_Albert_epilithic_algae_scraper,ylim=c(0,1),col="darkorange",pch=19,cex=1.5,
       xaxt="n",xlab="",ylab="Allele frequency",las=2,xlim=c(0,111000))
  points(frqZoom$POS,frqZoom$Frq_Albert_non_epilithic_algae_scraper,pch=8)
  points(frqZoom$POS[is.na(frqZoom$Frq_Albert_epilithic_algae_scraper)],rep(1,length=length(which(is.na(frqZoom$Frq_Albert_epilithic_algae_scraper)))),pch=0,col="darkorange")
  abline(v=refSubset$add)
  legend("topleft",legend=c("1 epilithic algae scraper","10 other cichlids  "),
         pch=c(19,8),col=c("darkorange","black"),bty="n",title="  Albert", title.adj=0)
  mtext(2,text = "H",las=2,line=3,at = 1,cex=1.1)
  
  # Kivu
  plot(frqZoom$POS,frqZoom$Frq_Kivu_epilithic_algae_scraper,ylim=c(0,1),col="darkorange",pch=19,cex=1.5,
       xaxt="n",xlab="",ylab="Allele frequency",las=2,xlim=c(0,111000))
  points(frqZoom$POS,frqZoom$Frq_Kivu_non_epilithic_algae_scraper,pch=8)
  points(frqZoom$POS[is.na(frqZoom$Frq_Kivu_epilithic_algae_scraper)],rep(1,length=length(which(is.na(frqZoom$Frq_Kivu_epilithic_algae_scraper)))),pch=0,col="darkorange")
  abline(v=refSubset$add)
  legend("topleft",legend=c("8 epilithic algae scraper","80 other cichlids  "),
         pch=c(19,8),col=c("darkorange","black"),bty="n",title="  Kivu", title.adj=0)
  mtext(2,text = "I",las=2,line=3,at = 1,cex=1.1)
  axis(1,at=refSubset$add+refSubset$length/2,refSubset$scaff,tick=F,line=1.3)
  axis(1,frqZoom$POS,labels=round(frqZoom$POS/1000000,digits=2))
  mtext(1,text = "Chromosome 11 position [Mb]",line=4)
  
  
}

dev.off()



