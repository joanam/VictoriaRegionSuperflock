#**************#
# Plot lake MDS 
#**************#

# Load libraries
require("SNPRelate")
require("Rtsne")
require("GWASTools")
require("readxl")

# Import GDS file for the analysis
gds <- snpgdsOpen(filename="D://Dropbox/victoriaGenomes/PCA/allGenomes.chr1-22.SNPs.minDP6.minGQ20.max0.5N.3masks.LDpruned.1kbThinned.gds")

# Import sample identifiers
id=read.gdsn(index.gdsn(gds,"sample.id"))


# Combine with sample information 
m<-as.data.frame(read_excel("D:/Dropbox/CichlidGenomesProject/planning/Cichlids_Sequenced_New.xlsx"))
info=m

# Separate Northern and Southern Generalists from the Victoria Radiation
m$Group[grepl(m$Species_Name_CleanDM, pattern="nubila swamp")]<-"LVRS Southern Generalists"
m$Group[grepl(m$Species_Name_CleanDM, pattern="ruby|flameback|cherry|lividus-like")]<-"LVRS Northern Generalists"


#### Prepare colours and symbols ####

m$Ecology[m$Genus=="Incertae sedis"]<-"?"
m$Ecology[m$sample=="109768"]<-"detritivore"
m$Genus[m$sample=="109768"]<-"Enterochromis"
m$Ecology[m$sample=="80805"]<-"?"
m$Genus[m$sample=="80805"]<-"Incertae sedis"
m$Ecology[m$sample%in%info$sample[info$Ecology2=="dwarfPredator"]]<-"dwarf predator"
m[m$sample=="130735","Ecology"]<-"?"
m[m$sample=="130735","Genus"]<-"Incertae sedis" # Probably Prognathochromis
m[m$Species=="demersal sp7"|m$Species=="demersal sp10"|m$Species=="piscivorous sp2","Group"]<-"Astatotilapia (Nile)"


# Extract species names
specname=m[match(id,m$sample),c("Species_Name_CleanDM")]
specname[m$sample=="109768"]<-"Enterochromis cinctus StE"


# Extract group memberships
groups=m[match(id,m$sample),c("Group")]
genGroup<-m[match(id,m$sample),c("genomicGroup")]
genus<-m[match(id,m$sample),"Genus"]

# simplify genera
genus[m[match(id,m$sample),"Lake"]!="Victoria" & m[match(id,m$sample),"Lake"]!="Kyoga"]<-"Incertae sedis"
genus[grepl(genus,pattern="Yssichromis")]<-"Yssichromis"
genus[grepl(genus,pattern="Pundamilia")]<-"Pundamilia"
genus[grepl(genus,pattern="Enterochromis")]<-"Enterochromis"
genus[grepl(genus,pattern="Paralabidochromis")]<-"Paralabidochromis"
genus[grepl(genus,pattern="Thoracochromis")]<-"Thoracochromis"
genus[grepl(genus,pattern="Throacochromis")]<-"Thoracochromis"
genus[grepl(genus,pattern="Gaurochromis")]<-"Gaurochromis"
genus[grepl(genus,pattern="Harpagochromis")]<-"Harpagochromis"
genus[grepl(genus,pattern="Psammochromis")]<-"Psammochromis"
genus[grepl(genus,pattern="Lipochromis")]<-"Lipochromis"
genus[grepl(genus,pattern="Neochromis")]<-"Neochromis"
genus[grepl(genus,pattern="Hoplotilapia")]<-"Hoplotilapia"
genus[grepl(genus,pattern="Ptyochromis")]<-"Ptyochromis"
genus[grepl(genus,pattern="Macropleurodus")]<-"Macropleurodus"
genus[grepl(genus,pattern="Harpago/Prognatho")]<-"Harpagochromis"

# Extract ecology
ecology=m[match(id,m$sample),c("Ecology")]
# Simplify ecology variable
ecology[grep("[?]",ecology)]="?"
ecology[grep(genGroup,pattern="doublestripe")]="dwarf predator"
ecology[grep(genGroup,pattern="dwarfPred")]="dwarf predator"
ecology[grep("might also eat fish",ecology)]="insectivore"
ecology[grep("oral crusher",ecology)]="oral snail crusher"
ecology[grep("probably NOT",ecology)]="?"
ecology[grep("oral sheller",ecology)]="oral snail sheller"
ecology[grep("oral sheller/detritivore",ecology)]="oral snail sheller"
ecology[grep("snail crusher",ecology)]="oral snail crusher"
ecology[grep("insect picker",ecology)]="insectivore"
ecology[grep("generalist-algae scraper",ecology)]="generalist"
ecology[grep("insectivore-pharyngeal-crusher",ecology)]="pharyngeal snail crusher"
ecology[grep("piscivore-plant-eater",ecology)]="piscivore"
ecology[grep("epiphytic algae scraper/insectivore",ecology)]="algae scraper/insectivore"
ecology[grep("detritivore/phytoplanktivore",ecology)]="detritivore"
ecology[grep("oral mollusc crusher",ecology)]="oral snail crusher"
ecology[grep("pharyngeal mollusc crusher",ecology)]="pharyngeal snail crusher"
ecology[grepl("Pundamilia",genus)&ecology=="zooplanktivore"]="reef zooplanktivore"
ecology[ecology=="zooplanktivore"]="pelagic zooplanktivore"

species=m[match(id,m$sample),c("Species")]


# Define colors for groups and symbols for ecologies
# Define color blind save colors from Wong et al. 2011 Nature Methods

# Colours
wongcol=apply(matrix(c(0,0,0,230,159,0,86,180,233,0,158,115,240,228,66,0,114,178,213,94,0,204,121,167),3)/255,2,function(x){rgb(x[1],x[2],x[3])})
colf=colorRampPalette(wongcol[c(1,6,3,4,5,2,7,8)])
colg=colf(length(unique(groups)))
groupcol=adjustcolor(colg[match(groups,unique(groups))],0.7)
# groupcol<-as.factor(groupcol)

cole=colf(length(unique(ecology)))
ecocol=adjustcolor(cole[match(ecology,unique(ecology))],0.7)

colgen=colf(length(unique(genus)))
gencol=adjustcolor(colgen[match(genus,unique(genus))],0.7)

# Symbols
pche=c(6,2,1,3,15,16,4,18,17,10,7,0,8,5,19,14,11,9)
ecopch=pche[match(ecology,unique(ecology))]
pchg=c(19,9,2,3,0,4,15,5,7,25,1,10,11,17,12,8,6,13,14,17,22,23)[1:length(unique(groups))]
grouppch=pchg[match(groups,unique(groups))]
pchgen=c(19,9,2,3,0,4,15,5,7,25,1,10,11,16,12,8,6,13,14,17,22,23,20,21,24)[1:length(unique(genus))]
genpch=pchgen[match(genus,unique(genus))]

#### Functions #####
plotMDS<-function(sub,poly=T,sub2="",letter=""){
  info=m[match(id,as.character(m$sample)),]
  
  infoSub<-info[sub,]
  
  # Compute distance matrix: IBS pairwise distances
  ibs=snpgdsIBS(gds,sample.id=id[sub],remove.monosnp = T,maf=0.05)
  
  # Run multidimensional scaling analysis on distance matrix
  mds=cmdscale(1-ibs$ibs,k=2)
  
  # Define coloring / points by groups / ecology
  cols=ecocol[sub]
  pchs=ecopch[sub]
  legends=c(unique(ecology[sub]))
  
  # Plot MDS
  plot(mds,xlab="MDS Factor 1",ylab="MDS Factor 2",col=cols,pch=pchs,bg="#00000066")
  
  # Draw convex hull among same ecology groups:
  ecoSub<-ecology[sub]
  if(poly){
    for(i in unique(ecoSub)[!("?"==unique(ecoSub))]){
      if(length(which(ecoSub==i))>1){
        ci=chull(mds[which(ecoSub==i&!infoSub$Group==sub2,),])
        polygon(mds[which(ecoSub==i&!infoSub$Group==sub2,)[ci],1],
                mds[which(ecoSub==i&!infoSub$Group==sub2,)[ci],2],
                col=NA,
                border=cols[match(i,ecoSub)])
      }
    }
  }
  if (poly&sub2!=""){
    for(i in unique(ecoSub)[!("?"==unique(ecoSub))]){
      if(length(which(ecoSub==i))>1){
        ci=chull(mds[which(infoSub$Group==sub2,),])
        polygon(mds[which(infoSub$Group==sub2,)[ci],1],
                mds[which(infoSub$Group==sub2,)[ci],2],
                col=NA,
                border="black")
      }
    }
  }
  mtext(side = 2,text = letter,las=2,line = 1.5,at = max(mds[,2])*1.03,cex=1.4)
  
}


#### Fig. 2 Manuscript figure ####

sub=(1:length(groups))[(groups=="LVRS Edward"| groups=="LVRS Albert" |
                          groups=="LVRS Saka" | groups=="LVRS Kagera"|
                          groups=="LVRS Kyoga" | groups=="LVRS Egypt" |
                          groups=="LVRS Sudan" | groups=="LVRS Turkana" |
                          groups=="LVRS Sahara") ]
sub=c(sub,sample((1:length(groups))[groups=="LVRS Victoria" ],size = 25))
sub=c(sub,sample((1:length(groups))[groups=="LVRS Kivu"],size = 25))
sub<-sort(sub)

# Compute distance matrix: IBS pairwise distances
ibs=snpgdsIBS(gds,sample.id=id[sub],remove.monosnp = T)

# Run multidimensional scaling analysis on distance matrix
mds=cmdscale(1-ibs$ibs,k=2)


# png(filename = "D:/Dropbox/victoriaGenomes/PCA/MDS_LVRS_Fig2.png",
#     width=200*(max(mds[,1])-min(mds[,1])),
#     height = 200*(max(mds[,2])-min(mds[,2])), units = 'cm', res = 300)

#{

# Define coloring / points by groups / ecology
info<-as.data.frame(cbind(groups,groupcol,grouppch))
info<-info[sub,]
info<-droplevels(info)
info$groupcol<-as.factor(info$groupcol)
levels(info$groupcol)<-c("#FF9F00","#FF6304","darkblue","red","#1DA180","#1DA180",
                         "#1DA180","#1DA180","deepskyblue","cornflowerblue","blue")

pchs<-as.integer(as.character(info$grouppch))

lakeCol<-factor(info$groupcol)
# levels(lakeCol)<-c("blue","slateblue","red","cornflowerblue","#FF6304","blue4","#FF9F00")
cols<-as.character(lakeCol)
legends<-unique(info)
legends<-rbind(legends[c(-3,-6),],legends[c(3,6),])

# Plot MDS
par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(4,4,1,1),xaxs="r",yaxs="r",mgp=c(1.5,0.5,0))
plot(-mds[,1],mds[,2],xlab="MDS Factor 1",ylab="MDS Factor 2",col=cols,pch=pchs)
legend("bottomright",
       substr(legends$groups,6,100),bty="n",pch=as.integer(as.vector(legends$grouppch)),
       col=as.character(legends$groupcol),
       ncol=2,cex=1)

# Draw convex hulls among lakes (excluding squamipinnis, Edward egg sneaker)
for(i in unique(groups[sub])[!("?"==unique(groups[sub]))]){
  if(length(which(groups[sub]==i))>1){
    mds_subset<-which(groups[sub]==i&!grepl(species[sub],pattern="squamipinnis|sp. Edward|egg sneaker|Mpanga"))
    ci=chull(mds[mds_subset,])
    polygon(-mds[mds_subset[ci],1],
            mds[mds_subset[ci],2],
            col=NA,
            border=cols[match(i,groups[sub])])
  }
}

# Label nubila
points(-mds[grepl(species[sub],pattern="nubila swamp"),1],
       mds[grepl(species[sub],pattern="nubila swamp"),2],pch=19,cex=0.5)
points(-mds[grepl(species[sub],pattern="cherry|green|lividus-like|flameback"),1],
       mds[grepl(species[sub],pattern="cherry|green|lividus-like|flameback"),2],pch=8,cex=0.5)
# }
# 
# dev.off()



### Fig. S4 MDS of different lakes ####

png(file="D:/Dropbox/victoriaGenomes/PCA/MDS_lakes_FigS6.png",width=850,height = 600)


par(mfrow=c(2,3),mar=c(4,4,1,1))

##### A LVRS incl populations further downstream of the Nile #### 

{
  sub=grep(groups,pattern = "LVRS")
  sub<-sub[-which(groups[sub]=="LVRS Rusizi")]
  
  # Compute PCA / MDS
  # Run principal components analysis
  # apca=snpgdsPCA(gds,sample.id=id[sub],remove.monosnp = T,)
  
  # Compute distance matrix: IBS pairwise distances
  ibs=snpgdsIBS(gds,sample.id=id[sub],remove.monosnp = T)
  
  # Run multidimensional scaling analysis on distance matrix
  mds=cmdscale(1-ibs$ibs,k=2)
  
  # Define coloring / points by groups / ecology
  info<-as.data.frame(cbind(groups,groupcol,grouppch))
  info<-info[sub,]
  info<-droplevels(info)
  info$groupcol<-as.factor(info$groupcol)
  levels(info$groupcol)<-c("#FF9F00","darkblue","#FF6304","red","#1DA180","#1DA180","#1DA180","deepskyblue","cornflowerblue","blue")
  cols<-as.character(info$groupcol)
  pchs<-as.integer(as.character(info$grouppch))
  
  legends<-unique(info)
  
  
  # Plot MDS
  plot(mds,xlab="MDS Factor 1",ylab="MDS Factor 2",col=cols,pch=pchs,
       bg="#00000066")
  legend("bottomright",
         substr(legends$groups,6,100),bty="n",pch=as.integer(as.vector(legends$grouppch)),col=as.character(legends$groupcol),
         pt.bg="#00000066",ncol=2)
  # Draw convex hull among same ecology groups:
  for(i in unique(groups[sub])[!("?"==unique(groups[sub]))]){
    if(length(which(groups[sub]==i))>1){
      ci=chull(mds[which(groups[sub]==i),])
      polygon(mds[which(groups[sub]==i)[ci],1],mds[which(groups[sub]==i)[ci],2],
              col=NA,
              border=cols[match(i,groups[sub])])
    }
  }
  mtext(side = 2,text = "A",las=2,line = 2.5,at = max(mds[,2])*1.03,cex=1.4)
}


##### B Western Lakes and Lower Nile #### 

{
  sub=(1:length(groups))[groups=="LVRS Edward"|groups=="LVRS Kivu"| groups=="LVRS Albert" | 
                           groups=="LVRS Saka" | groups=="LVRS Sudan" | groups =="LVRS Egypt" |
                           groups=="LVRS Boukou" | groups=="LVRS Turkana"]
  
  # Compute PCA / MDS
  # Run principal components analysis
  # apca=snpgdsPCA(gds,sample.id=id[sub],remove.monosnp = T,)
  
  # Compute distance matrix: IBS pairwise distances
  ibs=snpgdsIBS(gds,sample.id=id[sub],remove.monosnp = T)
  
  # Run multidimensional scaling analysis on distance matrix
  mds=cmdscale(1-ibs$ibs,k=2)
  
  # Define coloring / points by groups / ecology
  info<-as.data.frame(cbind(groups,groupcol,grouppch))
  info<-info[sub,]
  info<-droplevels(info)
  info$groupcol<-as.factor(info$groupcol)
  levels(info$groupcol)<-c("darkblue","#1DA180","#1DA180","#1DA180","deepskyblue","cornflowerblue","blue")
  cols<-as.character(info$groupcol)
  pchs<-as.integer(as.character(info$grouppch))
  
  legends<-unique(info)
  
  
  # Plot MDS
  plot(mds,xlab="MDS Factor 1",ylab="MDS Factor 2",col=cols,pch=pchs,
       bg="#00000066")
  legend("top",
         substr(legends$groups,6,100),bty="n",pch=as.integer(as.vector(legends$grouppch)),col=as.character(legends$groupcol),
         pt.bg="#00000066",ncol=2)
  # Draw convex hull among same ecology groups:
  for(i in unique(groups[sub])[!("?"==unique(groups[sub]))]){
    if(length(which(groups[sub]==i))>1){
      ci=chull(mds[which(groups[sub]==i),])
      polygon(mds[which(groups[sub]==i)[ci],1],mds[which(groups[sub]==i)[ci],2],
              col=NA,
              border=cols[match(i,groups[sub])])
    }
  }
  mtext(side = 2,text = "B",las=2,line = 3,at = max(mds[,2])*1.03,cex=1.4)
  
}


##### C Western Lakes #### 

{
  sub=(1:length(groups))[groups=="LVRS Edward"|groups=="LVRS Kivu"| groups=="LVRS Albert" | 
                           groups=="LVRS Saka"]
  
  # Compute PCA / MDS
  # Run principal components analysis
  # apca=snpgdsPCA(gds,sample.id=id[sub],remove.monosnp = T,)
  
  # Compute distance matrix: IBS pairwise distances
  ibs=snpgdsIBS(gds,sample.id=id[sub],remove.monosnp = T)
  
  # Run multidimensional scaling analysis on distance matrix
  mds=cmdscale(1-ibs$ibs,k=2)
  
  # Define coloring / points by groups / ecology
  info<-as.data.frame(cbind(groups,groupcol,grouppch))
  info<-info[sub,]
  info<-droplevels(info)
  info$groupcol<-as.factor(info$groupcol)
  levels(info$groupcol)<-c("darkblue","deepskyblue","cornflowerblue","blue")
  cols<-as.character(info$groupcol)
  pchs<-as.integer(as.character(info$grouppch))
  
  legends<-unique(info)
  
  
  # Plot MDS
  plot(mds,xlab="MDS Factor 1",ylab="MDS Factor 2",col=cols,pch=pchs,
       bg="#00000066")
  legend("bottomleft",
         substr(legends$groups,6,100),bty="n",pch=as.integer(as.vector(legends$grouppch)),col=as.character(legends$groupcol),
         pt.bg="#00000066",ncol=3)
  # Draw convex hull among same ecology groups:
  for(i in unique(groups[sub])[!("?"==unique(groups[sub]))]){
    if(length(which(groups[sub]==i))>1){
      ci=chull(mds[which(groups[sub]==i),])
      polygon(mds[which(groups[sub]==i)[ci],1],mds[which(groups[sub]==i)[ci],2],
              col=NA,
              border=cols[match(i,groups[sub])])
    }
  }
  mtext(side = 2,text = "C",las=2,line = 3,at = max(mds[,2])*1.03,cex=1.4)
  
}


##### D Edward, Albert and Kivu #### 

{
  sub=(1:length(groups))[(groups=="LVRS Edward"|groups=="LVRS Kivu"| groups=="LVRS Albert")]
  
  # Compute PCA / MDS
  # Run principal components analysis
  # apca=snpgdsPCA(gds,sample.id=id[sub],remove.monosnp = T,)
  
  # Compute distance matrix: IBS pairwise distances
  ibs=snpgdsIBS(gds,sample.id=id[sub],remove.monosnp = T)
  
  # Run multidimensional scaling analysis on distance matrix
  mds=cmdscale(1-ibs$ibs,k=2)
  
  # Define coloring / points by groups / ecology
  info<-as.data.frame(cbind(groups,groupcol,grouppch))
  info<-info[sub,]
  info$groupcol<-as.factor(info$groupcol)
  info$groupcol<-droplevels(info$groupcol)
  levels(info$groupcol)<-c("darkblue","cornflowerblue","blue")
  cols<-as.character(info$groupcol)
  pchs<-as.integer(as.character(info$grouppch))
  
  legends<-unique(info)
  
  
  # Plot MDS
  plot(mds,xlab="MDS Factor 1",ylab="MDS Factor 2",col=cols,pch=pchs,
       bg="#00000066")
  legend("bottomleft",
         substr(legends$groups,6,100),bty="n",pch=as.integer(as.vector(legends$grouppch)),col=as.character(legends$groupcol),
         pt.bg="#00000066")
  # Draw convex hull among same ecology groups:
  for(i in unique(groups[sub])[!("?"==unique(groups[sub]))]){
    if(length(which(groups[sub]==i))>1){
      ci=chull(mds[which(groups[sub]==i),])
      polygon(mds[which(groups[sub]==i)[ci],1],mds[which(groups[sub]==i)[ci],2],
              col=NA,
              border=cols[match(i,groups[sub])])
    }
  }
  mtext(side = 2,text = "D",las=2,line = 2.5,at = max(mds[,2])*1.03,cex=1.4)
  
}


##### E Victoria-Kyoga and Kagera #### 

{
  sub=(1:length(groups))[which(groups=="LVRS Victoria"|groups=="LVRS Kagera"|groups=="LVRS Kyoga"| groups=="LVRS Ikimba"| groups=="LVRS Burigi" | groups=="LVRS Nabugabo")]
  
  # Compute PCA / MDS
  # Run principal components analysis
  # apca=snpgdsPCA(gds,sample.id=id[sub],remove.monosnp = T,)
  
  # Compute distance matrix: IBS pairwise distances
  ibs=snpgdsIBS(gds,sample.id=id[sub],remove.monosnp = T)
  
  # Run multidimensional scaling analysis on distance matrix
  mds=cmdscale(1-ibs$ibs,k=2)
  
  # Define coloring / points by groups / ecology
  info<-as.data.frame(cbind(groups,groupcol,grouppch))
  info<-info[sub,]
  info$groupcol<-as.factor(info$groupcol)
  info$groupcol<-droplevels(info$groupcol)
  levels(info$groupcol)<-c("#FF9F00","#FF6304","red")
  cols<-as.character(info$groupcol)
  pchs<-as.integer(as.character(info$grouppch))
  
  legends<-unique(info)
  
  
  # Plot MDS
  plot(mds,xlab="MDS Factor 1",ylab="MDS Factor 2",col=cols,pch=pchs,
       bg="#00000066")
  legend("topright",
         substr(legends$groups,6,100),bty="n",pch=as.integer(as.vector(legends$grouppch)),col=as.character(legends$groupcol),
         pt.bg="#00000066")
  # Draw convex hull among same ecology groups:
  for(i in unique(groups[sub])[!("?"==unique(groups[sub]))]){
    if(length(which(groups[sub]==i))>1){
      ci=chull(mds[which(groups[sub]==i),])
      polygon(mds[which(groups[sub]==i)[ci],1],mds[which(groups[sub]==i)[ci],2],
              col=NA,
              border=cols[match(i,groups[sub])])
    }
  }
  mtext(side = 2,text = "E",las=2,line = 3,at = max(mds[,2])*1.03,cex=1.4)
  
}


##### F Victoria-Kyoga and Kagera but Victoria/Kyoga subsampled to match the size of Kagera ####

{
  sub=(1:length(groups))[which(groups=="LVRS Victoria"|groups=="LVRS Kyoga"| groups=="LVRS Ikimba"| groups=="LVRS Burigi" | groups=="LVRS Nabugabo")]
  # sub=c(sub[1:30],
  #       (1:length(groups))[groups=="LVRS Kagera"])
  sub=sort(c(sample(x = sub,size=length(which(groups=="LVRS Kagera"))),
             (1:length(groups))[groups=="LVRS Kagera"]))
  
  # Compute PCA / MDS
  # Run principal components analysis
  # apca=snpgdsPCA(gds,sample.id=id[sub],remove.monosnp = T,)
  
  # Compute distance matrix: IBS pairwise distances
  ibs=snpgdsIBS(gds,sample.id=id[sub],remove.monosnp = T)
  
  # Run multidimensional scaling analysis on distance matrix
  mds=cmdscale(1-ibs$ibs,k=2)
  
  # Define coloring / points by groups / ecology
  info<-as.data.frame(cbind(groups,groupcol,grouppch))
  info<-info[sub,]
  info$groupcol<-as.factor(info$groupcol)
  info$groupcol<-droplevels(info$groupcol)
  levels(info$groupcol)<-c("#FF9F00","#FF6304","red")
  cols<-as.character(info$groupcol)
  pchs<-as.integer(as.character(info$grouppch))
  
  legends<-unique(info)
  
  
  # Plot MDS
  plot(mds,xlab="MDS Factor 1",ylab="MDS Factor 2",col=cols,pch=pchs,
       bg="#00000066")
  legend("topright",
         substr(legends$groups,6,100),bty="n",pch=as.integer(as.vector(legends$grouppch)),col=as.character(legends$groupcol),
         pt.bg="#00000066")
  # Draw convex hull among same ecology groups:
  for(i in unique(groups[sub])[!("?"==unique(groups[sub]))]){
    if(length(which(groups[sub]==i))>1){
      ci=chull(mds[which(groups[sub]==i),])
      polygon(mds[which(groups[sub]==i)[ci],1],
              mds[which(groups[sub]==i)[ci],2],
              col=NA,
              border=cols[match(i,groups[sub])])
    }
  }
  mtext(side = 2,text = "F",las=2,line = 3,at = max(mds[,2])*1.03,cex=1.4)
  
}

dev.off() 




#### Fig S9: MDS of Victoria-Kyoga with two different symbol sets ####

png(file="D:/Dropbox/victoriaGenomes/PCA/MDS_Victoria-Kyoga.png",width=850,height = 500)

par(mfrow=c(1,2),mar=c(4,4,1,1),oma=c(0,0,0,0))

# Victoria-Kyoga Ecogroup
{
  sub=(1:length(groups))[groups=="LVRS Victoria"|groups=="LVRS Kyoga"]
  
  # Compute distance matrix: IBS pairwise distances
  ibs=snpgdsIBS(gds,sample.id=id[sub],remove.monosnp = T)
  
  # Run multidimensional scaling analysis on distance matrix
  mds=cmdscale(1-ibs$ibs,k=2)
  
  # Define coloring / points by groups / ecology
  info<-as.data.frame(cbind(ecology,ecocol,ecopch))
  info<-info[sub,]
  info<-droplevels(info)
  cols<-as.character(info$ecocol)
  pchs<-as.integer(as.character(info$ecopch))
  
  labels<-species[sub]
  
  legends<-unique(info)
  nubila<-grep(labels,pattern="swamp")
  northgen<-grep(labels,pattern="swamp")
  
  # Plot MDS
  plot(mds,xlab="MDS Factor 1",ylab="MDS Factor 2",col=cols,pch=pchs,
       bg="#00000066")
  for(i in unique(ecology[sub])[!("?"==unique(ecology[sub]))]){
    if(length(which(ecology[sub]==i))>1){
      ci=chull(mds[which(ecology[sub]==i),])
      polygon(mds[which(ecology[sub]==i)[ci],1],mds[which(ecology[sub]==i)[ci],2],
              col=NA,
              border=cols[match(i,ecology[sub])])
    }
  }
  legend("topright",col=as.character(legends$ecocol),bty="n",
         pch=as.integer(legends$ecopch),legend=legends$ecology,ncol=2,cex=0.8)
  
  # Show where nubila is and the Northern Generalists
  points(mds[nubila,],cex=0.5,pch=19)
  points(mds[northgen,],cex=0.5,pch=19)
  
  text(x = -0.026,y = 0.0342,labels = "A",las=2,xpd=T,cex=2)
}


# Victoria-Kyoga Genera
{
  # Define coloring / points by groups / ecology
  info<-as.data.frame(cbind(genus,gencol,genpch))
  info<-info[sub,]
  info<-droplevels(info)
  cols<-as.character(info$gencol)
  pchs<-as.integer(as.character(info$genpch))
  
  legends<-unique(info)
  
  
  # Plot MDS
  plot(mds,xlab="MDS Factor 1",ylab="MDS Factor 2",col=cols,pch=pchs,
       bg="#00000066")
  for(i in unique(genus[sub])[!("Incertae sedis"==unique(genus[sub]))]){
    if(length(which(genus[sub]==i))>1){
      ci=chull(mds[which(genus[sub]==i),])
      polygon(mds[which(genus[sub]==i)[ci],1],mds[which(genus[sub]==i)[ci],2],
              col=NA,
              border=cols[match(i,genus[sub])])
    }
  }
  legend("topright",col=as.character(legends$gencol),pch=as.integer(legends$genpch),
         legend=legends$genus,ncol=2,cex=0.8,bty="n")
  text(x = -0.026,y = 0.0342,labels = "B",las=2,xpd=T,cex=2)
  
}

dev.off()




#### Subsampled to 4 Victoria/Kyoga, plus all Northern, all Southern Generalists and 4 Edward inds ####

LVseparateCount<-0
LVcentral<-0
both<-0

pdf(file="D:/Dropbox/victoriaGenomes/PCA/MDS_subsampled_Generalists_Edward_Kivu.pdf",
    height=7,width = 7)
for(i in 1:1000)
{
  # Downsample to four individuals
  subVic=which(groups=="LVRS Victoria"|groups=="LVRS Kyoga")
  subKiv=which(groups=="LVRS Kivu")
  
  subSouthGen=which(groups=="LVRS Southern Generalists")
  sub=sort(c(sample(x = subVic,size=4),sample(x = subKiv,size=4),
             sample(x=subSouthGen,size=4),
             grep(groups,pattern="Northern Generalists")))
  
  # Compute distance matrix: IBS pairwise distances
  ibs=snpgdsIBS(gds,sample.id=id[sub],remove.monosnp = T)
  
  # Run multidimensional scaling analysis on distance matrix
  mds=cmdscale(1-ibs$ibs,k=2)
  
  # Define coloring / points by groups / ecology
  info<-as.data.frame(cbind(groups,groupcol,grouppch))
  info<-info[sub,]
  info[info$groups=="LVRS Kyoga",]<-info[info$groups=="LVRS Victoria",][1,]
  info$groupcol<-as.factor(info$groupcol)
  info$groupcol<-droplevels(info$groupcol)
  
  legends<-unique(info)
  
  info$MDS1<-mds[,1]
  info$MDS2<-mds[,2]
  
  
  # levels(info$groupcol)<-c("#FF9F00","#FF6304","red")
  cols<-as.character(info$groupcol)
  pchs<-as.integer(as.character(info$grouppch))
  
  
  # Check if Victoria is outside all other groups:
  require(sf)
  vi=chull(mds[info$groups=="LVRS Victoria",])
  vi=c(vi,vi[1])
  viP<-st_polygon(list(mds[which(info$groups=="LVRS Victoria")[vi],1:2]))
  ng=chull(mds[info$groups=="LVRS Northern Generalists",])
  ng=c(ng,ng[1])
  ngP<-st_polygon(list(mds[which(info$groups=="LVRS Northern Generalists")[ng],1:2]))
  
  minX<-min(info[info$groups=="LVRS Victoria" | info$groups=="LVRS Northern Generalists","MDS1"])
  maxX<-max(info[info$groups=="LVRS Victoria" | info$groups=="LVRS Northern Generalists","MDS1"])
  minY<-min(info[info$groups=="LVRS Victoria" | info$groups=="LVRS Northern Generalists","MDS2"])
  maxY<-max(info[info$groups=="LVRS Victoria" | info$groups=="LVRS Northern Generalists","MDS2"])
  
  # plot(viP,col="orange",xlim=c(minX,maxX),ylim=c(minY,maxY))
  # plot(ngP,col="yellow",add=T)
  
  # Compute overlap
  overlap<-st_intersection(viP,ngP)
  if(length(overlap)==0){
    overlap = "VR separate"
    LVseparateCount=LVseparateCount+1
  }else{
    overlap = "VR overlapping"
  }
  
  # Is Victoria Radiation more in the middle than Northern Genearlists in the first two axes?
  mds1<-FALSE
  if(mean(info[info$groups=="LVRS Victoria","MDS1"])<0 & 
     mean(info[info$groups=="LVRS Victoria","MDS1"])>
     mean(info[info$groups=="LVRS Northern Generalists","MDS1"])){
    mds1<-TRUE
  }else if(mean(info[info$groups=="LVRS Victoria","MDS1"])>0 & 
           mean(info[info$groups=="LVRS Victoria","MDS1"])<
           mean(info[info$groups=="LVRS Northern Generalists","MDS1"])){
    mds1<-TRUE
  }
  mds2<-FALSE
  if(mean(info[info$groups=="LVRS Victoria","MDS2"])<0 & 
     mean(info[info$groups=="LVRS Victoria","MDS2"])>
     mean(info[info$groups=="LVRS Northern Generalists","MDS2"])){
    mds2<-TRUE
  }else if(mean(info[info$groups=="LVRS Victoria","MDS2"])>0 & 
           mean(info[info$groups=="LVRS Victoria","MDS2"])<
           mean(info[info$groups=="LVRS Northern Generalists","MDS2"])){
    mds2<-TRUE
  }  
  
  LVpos<-""
  if(mds1 & mds2){
    LVpos<-"central"
    LVcentral=LVcentral+1
  }else{
    LVpos<-"outside"
  }
  
  # Plot MDS
  plot(mds[,1:2],xlab="MDS Factor 1",ylab="MDS Factor 2",col=cols,pch=pchs,
       bg="#00000066")
  legend("topright",
         substr(legends$groups,6,100),bty="n",
         pch=as.integer(as.vector(legends$grouppch)),
         col=as.character(legends$groupcol),
         pt.bg="#00000066")
  mtext(side=3,text=paste(overlap,LVpos))
  # Draw convex hull among same ecology groups:
  # Get represented groups:
  for(i in unique(info$groupcol)){
    ci=chull(mds[info$groupcol==i,])
    polygon(mds[which(info$groupcol==i)[ci],1],
            mds[which(info$groupcol==i)[ci],2],
            col=NA,border=cols[match(i,info$groupcol)])
  }
  
  if(LVpos=="central" & overlap=="VR separate")
    both=both+1
  
}
dev.off()


