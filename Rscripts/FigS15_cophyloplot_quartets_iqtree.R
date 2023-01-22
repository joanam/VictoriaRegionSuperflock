require("phytools")
require("ape")

# Read in the two trees you want to compare
quartets<-read.tree("D:/Dropbox/victoriaGenomes/phylogenetics/jo_quartet_concordance.tree")
iqtree<-read.tree("D:/Dropbox/victoriaGenomes/Dstats/forMatt/LVRS.one.sample.per.group.nwk")

# Correct species labels
iqtree$tip.label[166]<-"Astatotilapia_burtoni"  # instead of Outgroup
iqtree$tip.label[iqtree$tip.label=="Xystichromis_sp_flameback"]<-"Astatotilapia_sp_flameback"
iqtree$tip.label[iqtree$tip.label=="Xystichromis_sp_ruby_green"]<-"Haplochromis_sp_rubygreen"
iqtree$tip.label[iqtree$tip.label=="Incertae_sedis_sp_cherry_fin"]<-"Astatotilapia_sp_cherry_fin"
quartets$tip.label[quartets$tip.label=="Xystichromis_sp_flameback"]<-"Astatotilapia_sp_flameback"
quartets$tip.label[quartets$tip.label=="Xystichromis_sp_ruby_green"]<-"Haplochromis_sp_rubygreen"
quartets$tip.label[quartets$tip.label=="Incertae_sedis_sp_cherry_fin"]<-"Astatotilapia_sp_cherry_fin"

# ladderise both trees
iqtree<-ladderize(iqtree)
quartets<-ladderize(quartets)

# Add branch lengths to quartets tree
quartets$edge.length<-rep(0.01,length(quartets$edge.length))

# Get a dataframe with the association (here just twice the same individual)
assoc<-cbind(quartets$tip.label, quartets$tip.label)

# Get the cophyloplot object
obj<-cophylo(iqtree,quartets,assoc=assoc,cex=0.1)

# Genus colour:
genusCol<-data.frame("species"=assoc[,1])
require(tidyverse)
genusCol$genus<-str_split_fixed(genusCol$species,pattern="_",n=2)[,1]
genusCol$genus[c(1:36,111:115,166)]<-NA
genusCol$genus[genusCol$genus=="Incertae"]<-NA
genusCol$genus[genusCol$genus=="Burigi"]<-NA
genusCol$genus[genusCol$genus=="Ikimba"]<-NA
genusCol$genus[genusCol$genus=="Double"]<-"double-stripe group"
genusCol$genus<-as.factor(genusCol$genus)
colFanTree<-as.data.frame(cbind("genus"=levels(factor(sortedSampleInfo$genusCol)),"color"=palGenera))
colFanTree[22,1]<-"Tridontochromis" # remove "?"
colFanTree[15,1]<-"Platytaeniodus" # remove "?"
genusCol<-merge(genusCol,colFanTree,by="genus",all.x = T,sort = F)
rownames(genusCol)<-genusCol$species
genusCol<-genusCol[assoc[,1],]
genusCol[is.na(genusCol$color),"color"]<-"lightgrey"

par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(0,0,0,0))

# Make a first cophyloplot
plot.cophylo(obj,link.type="curved",link.lwd=3,fsize=0.3,pts=F,ftype="reg",
             link.lty="solid",link.col=make.transparent(genusCol$color,0.9))

# Add legends
legend("bottomleft",legend=colFanTree$genus,
       fill=colFanTree$color,bty="n",title="Genus",cex=0.4,title.adj = 0)
