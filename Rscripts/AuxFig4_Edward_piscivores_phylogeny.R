library(ape)
library(phyloch)

# Read in the phylogeny with H. kimondo
tree<-read.tree("D:/Dropbox/victoriaGenomes/phylogenetics/edw_kivu_vic.chr1-22.minDP3.max0.5N.masks.thin1kb.phylip.treefile")

require("readxl")
sampleInfo<-as.data.frame(read_excel("D:/Dropbox/Manuscripts/LVRSgenomes/secondReviewRound/Auxiliary/AuxTable1_Sample_List.xlsx"))

outgroups<-c("81342","81343","81344","81345")
tree<-root(tree,outgroup = outgroups,edgelabel = T)
tree<-drop.tip(tree,c(outgroups,"81346","81027"))

#Plot the tree
par(mfrow=c(1,1),mar=c(0.5,0.5,0.5,0.5))

plotTree <- function(tree,title="",showLab=T,nodeLab=T,xlimm=0.05,tNode=70,size=0.5) {
  treetips<-as.character(tree$tip.label)
  treetipsorder=cbind(treetips,"treenames"=treetips,"indexx"=c(100:(length(treetips)+99)))
  mergedSampleInfo<-merge(sampleInfo,treetipsorder,by.x="sample",by.y="treetips",all.y=T)
  sortedSampleInfo<-mergedSampleInfo[order(mergedSampleInfo$indexx),]
  sortedSampleInfo<-droplevels(sortedSampleInfo)
  sortedSampleInfo<-unique(x=sortedSampleInfo)
  sortedSampleInfo<-cbind(sortedSampleInfo,"group_color"=sortedSampleInfo$Lake)
  sortedSampleInfo$group_color<-as.factor(sortedSampleInfo$group_color)
  colors<-colorRampPalette(c("black","red","grey","blue","orange"))
  levels(sortedSampleInfo$group_color)[levels(sortedSampleInfo$group_color)=="Edward"]<-"blue"
  levels(sortedSampleInfo$group_color)[levels(sortedSampleInfo$group_color)=="Kivu"]<-"darkblue"
  levels(sortedSampleInfo$group_color)[levels(sortedSampleInfo$group_color)=="Victoria"]<-"darkorange"
  levels(sortedSampleInfo$group_color)[levels(sortedSampleInfo$group_color)=="Kyoga"]<-"darkorange"
  levels(sortedSampleInfo$group_color)[levels(sortedSampleInfo$group_color)==""]<-"black"
  
  plot(tree,show.tip=F,x.lim=xlimm,main=title,edge.width = 2,cex=size,
       tip.color = as.character(sortedSampleInfo$group_color))
  tiplabels(pch = 19,cex=size,col=as.character(sortedSampleInfo$group_color))
  if(showLab) tiplabels(adj=-0.05,frame="none",text=paste(sortedSampleInfo$sample,
                                                          sortedSampleInfo$Genus,
                                                          sortedSampleInfo$Species),
                        cex=size,col=as.character(sortedSampleInfo$group_color))
  if(nodeLab) nodelabels(tree$node.label,cex=ifelse(as.integer(tree$node.label)>tNode,0.5,0.00001),frame="n",bg="white",adj=c(1,-0.5))
}

# Plot the phylogeny with H. kimondo
plotTree(tree,xlimm=0.005,showLab=T,nodeLab=T,tNode=90,size=0.4)

