# Mitochondrial versus nuclear phylogeny
require("ape")
require("phyloch")
require("phytools")

# Read in sample information
require("readxl")
sampleInfo<-as.data.frame(read_excel("D:/Dropbox/victoriaGenomes/Cichlids_Sequenced_New.xlsx"))
rownames(sampleInfo)<-sampleInfo$sample
  
# Read in tree files:
mtTree<-read.tree("D:/Dropbox/victoriaGenomes/mtDNA/genomes2020.mtDNA.max0.5N.minDP5.treefile")
nucTree<-read.tree("D:/Dropbox/victoriaGenomes/phylogenetics/allGenomes.1kbThinned.withSRA.chr1-22.alluaOut.noSquami.phylip.varsites.phy.treefile")

# Drop taxa missing in one of the trees:
nucTree<-drop.tip(nucTree,tip = nucTree$tip.label[!nucTree$tip.label%in%mtTree$tip.label])
mtTree<-drop.tip(mtTree,tip = mtTree$tip.label[!mtTree$tip.label%in%nucTree$tip.label])

# Remove alluaudi and egg sneaker samples
alluaudi<-c("103218","106754","10685","10719","131344","131201","105845","103881","81711")
nucTree<-drop.tip(nucTree,tip = alluaudi)
mtTree<-drop.tip(mtTree,tip = alluaudi)

# Root the trees on burtoni, then remove the outgroup
nucTree<-root(nucTree,outgroup = c("Aburtoni","131282"),edgelabel = T)
mtTree<-root(mtTree,outgroup = c("Aburtoni","131282"),edgelabel = T)

nucTree<-drop.tip(nucTree,tip = c("Aburtoni","131282"))
mtTree<-drop.tip(mtTree,tip = c("Aburtoni","131282"))

# Ladderize both trees
nucTree<-ladderize(nucTree)
mtTree<-ladderize(mtTree)

# Rename tip labels
nucTree$tip.label<-paste(nucTree$tip.label,sampleInfo[nucTree$tip.label,"Species_Name_CleanDM"])
mtTree$tip.label<-paste(mtTree$tip.label,sampleInfo[mtTree$tip.label,"Species_Name_CleanDM"])

rownames(sampleInfo)<-paste(sampleInfo$sample,sampleInfo$Species_Name_CleanDM)


# Get a dataframe with the association (here just twice the same individual)
assoc<-cbind(nucTree$tip.label, nucTree$tip.label)

# Order the sample info dataframe
sampleInfo<-sampleInfo[nucTree$tip.label,]

# Get the cophyloplot object
obj<-cophylo(nucTree,mtTree,assoc=assoc,cex=0.1,rotate=F)

# Lake colour:
lakeCol<-data.frame("lake"=assoc[,1])

require(tidyverse)

lakeCol$col<-as.factor(sampleInfo$Lake)
levels(lakeCol$col)<-c("cornflowerblue","#1DA180","orange","darkred","black","blue","#1DA180",
                       "orange","red","darkblue","darkblue","orange","black","slateblue",
                       "orange","black","slateblue","#1DA180","#1DA180","orange")
lakeCol$col<-as.character(lakeCol$col)
lakeCol$col[sampleInfo$Group=="Astatotilapia (Nile)"&!is.na(sampleInfo$Group)]<-"#000047"
lakeCol$col[grepl(sampleInfo$Species_Name_CleanDM,pattern="pharyngalis")]<-"#00005A"
lakeCol$col[sampleInfo$sample%in%c("130716","130715","130734","130707")]<-"darkorange"
lakeCol$col[grepl(sampleInfo$Species_Name_CleanDM,pattern="nubila swamp")]<-"#F4C430"
rownames(lakeCol)<-paste(sampleInfo$sample,sampleInfo$Species_Name_CleanDM)

# Lake colours tree order
nucTips<-obj$trees[[1]]$tip.label
mtDNATips<-obj$trees[[2]]$tip.label
nucLakes<-lakeCol[nucTips,]
mtLakes<-lakeCol[mtDNATips,]


# Make a cophyloplot for Aux. Fig. 1
pdf("D:/Dropbox/victoriaGenomes/phylogenetics/Cophyloplot_mtDNA_nuclear.pdf",width = 10,height=25)

  par(mfrow=c(1,1),mar=c(1,1,1,1),xpd=T)

  plot.cophylo(obj,link.type="curved",link.lwd=4,fsize=0.3,pts=F,ftype="reg",
               link.lty="solid",link.col=make.transparent(lakeCol$col,0.9))
  
  
  # Add bootstrap support symbols to both trees
  nodelabels.cophylo(pch=15,col=ifelse(as.integer(obj$trees[[1]]$node.label)==100,"grey",
                                ifelse(as.integer(obj$trees[[1]]$node.label)>90,"grey",NA)),
                                cex=ifelse(as.integer(obj$trees[[1]]$node.label)==100,0.7,
                                       ifelse(as.integer(obj$trees[[1]]$node.label)>90,0.4,NA)))
  nodelabels.cophylo(pch=15,col=ifelse(as.integer(obj$trees[[2]]$node.label)==100,"grey",
                                ifelse(as.integer(obj$trees[[2]]$node.label)>90,"grey",NA)),
                     cex=ifelse(as.integer(obj$trees[[2]]$node.label)==100,0.7,
                                       ifelse(as.integer(obj$trees[[2]]$node.label)>90,0.4,NA)),
                     which="right")
  
  # Add tip labels with lake colours
  tiplabels.cophylo(pch=19,cex=0.5,col=as.character(nucLakes$col))
  tiplabels.cophylo(pch=19,cex=0.5,col=as.character(mtLakes$col),which = "right")
  
  legend(x=0.35, y=0.6, pch=19,bty="n",title="Colours",title.adj=0.1,cex=0.6,
         col=c("darkred","#000047","black","#1DA180","cornflowerblue",
               "blue","darkblue","slateblue","#F4C430","darkorange","orange","red"),
         legend=c("Congolese lineage","Upper Nile lineage","Eastern","Pliocene Nile","Albert",
                "Edward","Kivu","Saka","A. nubila","Northern Generalists","Victoria","Kagera"))
  
  legend(x=0.35, y=0.45,pch=15,bty="n",title.adj=0.1,cex=0.6,
         legend=c("100",">90"),title="Bootstrap support",col="grey")
  
dev.off()


