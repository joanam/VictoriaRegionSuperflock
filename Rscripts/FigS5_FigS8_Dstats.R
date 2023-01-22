### D stats manuscript figures ####

#### Fig. S5A Upper Nile ancestry  ####

dstats<-read.table("D:/Dropbox/victoriaGenomes/Dstats/Dstats_UpperNile_Congo.results",header=T,sep="\t")
dstats$low<-dstats$D-3*(dstats$D/dstats$z)
dstats$high<-dstats$D+3*(dstats$D/dstats$z)

dstats$P1[dstats$P1=="nubila_swamp"]<-"Southern Generalists"
dstats$P1[dstats$P1=="Northern_generalists"]<-"Northern Generalists"

dstats<-dstats[order(dstats$D),]

par(xaxs="i",yaxs="i",mar=c(9,4,1,1))
plot(1:12,dstats$D[dstats$P3=="Edward-Nile"],ylim=c(0,0.45),
     xlim=c(0,13),xlab="",xaxt="n",ylab="D statistic",
     col=c("orange","orange","orange","red","#1DA180","#1DA180","slateblue","cornflowerblue","blue","#1DA180","darkblue","#1DA180"),pch=19,cex=1)
arrows(x0=1:12, y0=dstats$low[dstats$P3=="Edward-Nile"], 
       x1=1:12, y1=dstats$high[dstats$P3=="Edward-Nile"],lwd=2,code=0)
arrows(x0=1:12+0.2, y0=dstats$low[dstats$P3=="Kivu-Nile"], 
       x1=1:12+0.2, y1=dstats$high[dstats$P3=="Kivu-Nile"],lwd=2,code = 0)
points(1:12,dstats$D[dstats$P3=="Edward-Nile"],ylim=c(0,0.5),
       col=c("orange","orange","orange","red","#1DA180","#1DA180","slateblue","cornflowerblue","blue","#1DA180","darkblue","#1DA180"),pch=19,cex=1)
points(1:12+0.2,dstats$D[dstats$P3=="Kivu-Nile"],ylim=c(0,0.3),pch=17,cex=1,
       col=c("orange","orange","orange","red","#1DA180","#1DA180","slateblue","cornflowerblue","blue","#1DA180","darkblue","#1DA180"))
axis(at = 1:12,side = 1,labels = dstats$P1[dstats$P3=="Edward-Nile"],las=2)
abline(h=0)
legend("topleft",bty="n",pch=c(19,17),
       legend=c(as.expression(bquote(italic('T. pharyngalis'))),as.expression(bquote(italic('H. gracilior')))),title="P3",title.adj = 0.08)

#### Fig. S5B: Upper Nile and Congolese ancestry blocks ####

# size distribution using BBAA (LVRS like Congo) and BABA (LVRS like Nile) patterns

setwd("D:/Dropbox/victoriaGenomes/Dstats/")
# Window sizes
w<-10000
j=0.2

# Species
sp<-c("nubilaSwamp","Victoria","Kagera","Kivu","Edward","Albert","Saka","Turkana","Sudan","Egypt","Boukou")

# Species
sp<-c("nubilaSwamp","Victoria","Kivu","Edward","Turkana","Sudan","Egypt")

congo<-data.frame(matrix(ncol=6,nrow=length(sp)))
row.names(congo)<-sp
nile<-data.frame(matrix(ncol=6,nrow=length(sp)))
row.names(nile)<-sp

legend.labels<-sp
legend.labels[1]<-"A. nubila"

for(i in sp){
  # Read input data
  dat<-read.csv(paste0("BABA_BBAA.",i,"-Con-Nil-allu.10kb.out"),
                colClasses = c("character","integer","numeric")[c(1,2,2,2,2,2,2,3,3,3,3,3,3,3)],na.strings = c("NA","NA2"))
  
  # Calculate BABA-proportion
  d<-data.frame(dat[,c(1:4)],prop=dat$BABA/(dat$BABA+dat$BBAA))
  
  # Assign proportions to states (0:mixed, 1:Upper Nile, 2: Congolese, 3: Missing)
  f<-data.frame(d,state=c(0,1,2)[1+1*(d$prop>1-j)+2*(d$prop<j)])
  f$state[is.na(f$state)]<-3
  
  # Remove single empty blocks
  f<-f[unlist(mapply(FUN=rep.int,list(apply(cbind(rle(f$state)$lengths==1,rle(f$state)$values==3),1,sum)!=2),times=list(rle(f$state)$lengths))),]
  
  # Indices for start and end of contingous blocks
  # Upper Nile
  unstart<-c(1,which(f$state[-1]-f$state[-length(f$state)]!=0)+1)[rle(f$state)$values==1]
  unend<-which(f$state[-1]-f$state[-length(f$state)]!=0)[rle(f$state)$values==1]
  
  # Congolese
  cstart<-c(1,which(f$state[-1]-f$state[-length(f$state)]!=0)+1)[rle(f$state)$values==2]
  cend<-which(f$state[-1]-f$state[-length(f$state)]!=0)[rle(f$state)$values==2]
  
  # Calculate block sizes
  unblocksizes<-f$end[unend]-f$start[unstart]
  unblocksizes<-unblocksizes[unblocksizes>0]
  if(any(is.na(unblocksizes))){unblocksizes<-unblocksizes[-which(is.na(unblocksizes))]}
  
  cblocksizes<-f$end[cend]-f$start[cstart]
  cblocksizes<-cblocksizes[cblocksizes>0]
  if(any(is.na(cblocksizes))){cblocksizes<-cblocksizes[-which(is.na(cblocksizes))]}
  
  # Get barplot values
  c<-hist(cblocksizes,breaks=c(seq(w/2,50000+w/2,by=w),max(c(unblocksizes,cblocksizes))+w/2),plot=F)
  n<-hist(unblocksizes,breaks=c(seq(w/2,50000+w/2,by=w),max(c(unblocksizes,cblocksizes))+w/2),plot=F)
  
  congo[i,]<-c$counts/sum(c$counts)
  nile[i,]<-n$counts/sum(n$counts)
}

# Plot Upper Nile ancestry block size distribution for all species together
b<-barplot(as.matrix(nile),beside = T,xaxt="n",ylim=c(0,1),
           border=F,col=c("orange","darkorange","darkblue","blue","#1DA180","#1DA180","#1DA180"),
           xlab=paste("Upper Nile block size [kb]",sep=""),
           ylab="proportion of Upper Nile ancestry blocks",main="")
axis(1,b[-3+7*(1:(length(b)/7))],c("0-10","10-20","20-30","30-40","40-50",">60"))
legend("topright",fill=c("orange","darkorange","darkblue","blue","#1DA180","#1DA180","#1DA180"),
       legend=legend.labels,cex=1,bty="n")

# Plot Congolese ancestry block distribution for all species together
b<-barplot(as.matrix(congo),beside = T,xaxt="n",ylim=c(0,1),
           border=F,col=c("orange","darkorange","darkblue","blue","#1DA180","#1DA180","#1DA180"),
           xlab=paste("Congolese ancestry block size [kb]",sep=""),
           ylab="proportion of Congolese ancestry blocks",main="")
axis(1,b[-3+7*(1:(length(b)/7))],c("0-10","10-20","20-30","30-40","40-50",">60"))



#### Fig. S5D-I: Lake Correlations ############################

# Correlation of BABA/(BBAA+BABA) proportion in 10 kb windows

setwd("D:/Dropbox/victoriaGenomes/Dstats/")
j=0.2

getAncProp<-function(lake="lake"){
  # read in dataset
  vicD<-read.csv(paste0("BABA_BBAA.",lake,"-Con-Nil-allu.10kb.out"),
                 colClasses = c("character","integer","numeric")[c(1,2,2,2,2,2,2,3,3,3,3,3,3,3)],na.strings = c("NA","NA2"))
  vicD<-vicD[!is.na(vicD$ABBA),]
  
  # filter out windows with too low number of usable SNPs
  vicD<-vicD[vicD$ABBA+vicD$BABA+vicD$BBAA>5,]
  vicD<-data.frame(vicD[,c(1:4)],prop=vicD$BABA/(vicD$BABA+vicD$BBAA))
  
  # Assign proportions to states (0:mixed, 1:Upper Nile, 2: Congolese, 3: Missing)
  vicD<-data.frame(vicD,state=c(0,1,2)[1+1*(vicD$prop>1-j)+2*(vicD$prop<j)])
  vicD$state[is.na(vicD$state)]<-3
  
  vicD<-vicD[!is.na("scaffold"),]
  return(vicD)
}

# Get different datasets
vicD<-getAncProp("Victoria")
edwD<-getAncProp("Edward")
albD<-getAncProp("Albert")
turD<-getAncProp("Turkana")
sudD<-getAncProp("Sudan")
nubD<-getAncProp("nubilaSwamp")

# Combine datasets pairwise
datVicEdw<-merge(vicD[c("scaffold","start","prop","state")],
           edwD[c("scaffold","start","prop","state")],by=c("scaffold","start"))
names(datVicEdw)[3:6]<-c("nileProp_Victoria","state_Victoria","nileProp_Edward","state_Edward")

datVicNub<-merge(vicD[c("scaffold","start","prop","state")],
           nubD[c("scaffold","start","prop","state")],by=c("scaffold","start"))
names(datVicNub)[3:6]<-c("nileProp_Victoria","state_Victoria",
                   "nileProp_nubila","state_nubila")

datVicTur<-merge(vicD[c("scaffold","start","prop","state")],
           turD[c("scaffold","start","prop","state")],by=c("scaffold","start"))
names(datVicTur)[3:6]<-c("nileProp_Victoria","state_Victoria","nileProp_Turkana","state_Turkana")

datEdwTur<-merge(edwD[c("scaffold","start","prop","state")],
           turD[c("scaffold","start","prop","state")],by=c("scaffold","start"))
names(datEdwTur)[3:6]<-c("nileProp_Edward","state_Edward","nileProp_Turkana","state_Turkana")

datVicSud<-merge(vicD[c("scaffold","start","prop","state")],
           sudD[c("scaffold","start","prop","state")],by=c("scaffold","start"))
names(datVicSud)[3:6]<-c("nileProp_Victoria","state_Victoria","nileProp_Sudan","state_Sudan")

datTurSud<-merge(turD[c("scaffold","start","prop","state")],
           sudD[c("scaffold","start","prop","state")],by=c("scaffold","start"))
names(datTurSud)[3:6]<-c("nileProp_Turkana","state_Turkana","nileProp_Sudan","state_Sudan")

datEdwAlb<-merge(edwD[c("scaffold","start","prop","state")],
                 albD[c("scaffold","start","prop","state")],by=c("scaffold","start"))
names(datEdwAlb)[3:6]<-c("nileProp_Edward","state_Edward","nileProp_Albert","state_Albert")


# Make contour plots

require(ggplot2)
require(gridExtra)

g1<-ggplot(datVicNub, aes(x = nileProp_Victoria, y = nileProp_nubila)) +
  geom_point(size=0.2) +
  geom_density_2d_filled(alpha = 0.8) +
  theme_classic() + xlab("Victoria Radiation") + ylab("A. nubila") +
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold"),
        plot.tag.position = c(0.01,1)) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  geom_abline(
    mapping = NULL,
    data = NULL,
    intercept = 0, slope=1,
    na.rm = FALSE,
    show.legend = NA, colour="white"
  ) + labs(tag = "D") 
  


g2<-ggplot(datVicEdw, aes(x = nileProp_Victoria, y = nileProp_Edward)) +
  geom_point(size=0.2) +
  geom_density_2d_filled(alpha = 0.8) +
  theme_classic() +xlab("Victoria Radiation") + ylab("Lake Edward cichlids") +
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold"),
        plot.tag.position = c(0.01,1)) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  geom_abline(
    mapping = NULL,
    data = NULL,
    intercept = 0, slope=1,
    na.rm = FALSE,
    show.legend = NA, colour="white"
  )+ labs(tag = "E") 

g3<-ggplot(datVicTur, aes(x = nileProp_Victoria, y = nileProp_Turkana)) +
  geom_point(size=0.2) +
  geom_density_2d_filled(alpha = 0.8) +
  theme_classic() +xlab("Victoria Radiation") + ylab("Lake Turkana cichlids") +
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold"),
        plot.tag.position = c(0.01,1)) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  geom_abline(
    mapping = NULL,
    data = NULL,
    intercept = 0, slope=1,
    na.rm = FALSE,
    show.legend = NA, colour="white"
  )+ labs(tag = "F") 

g4<-ggplot(datEdwAlb, aes(x = nileProp_Edward, y = nileProp_Albert)) +
  geom_point(size=0.2) +
  geom_density_2d_filled(alpha = 0.8) +
  theme_classic() +xlab("Lake Edward cichlids") + ylab("Lake Albert cichlids") +
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold"),
        plot.tag.position = c(0.01,1)) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  geom_abline(
    mapping = NULL,
    data = NULL,
    intercept = 0, slope=1,
    na.rm = FALSE,
    show.legend = NA, colour="white"
  )+ labs(tag = "G") 

g5<-ggplot(datEdwTur, aes(x = nileProp_Edward, y = nileProp_Turkana)) +
  geom_point(size=0.2) +
  geom_density_2d_filled(alpha = 0.8) +
  theme_classic() +xlab("Lake Edward cichlids") + ylab("Lake Turkana cichlids") +
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold"),
        plot.tag.position = c(0.01,1)) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  geom_abline(
    mapping = NULL,
    data = NULL,
    intercept = 0, slope=1,
    na.rm = FALSE,
    show.legend = NA, colour="white"
  )+ labs(tag = "H") 

g6<-ggplot(datTurSud, aes(x = nileProp_Turkana, y = nileProp_Sudan)) +
  geom_point(size=0.2) +
  geom_density_2d_filled(alpha = 0.8) + 
  theme_classic() +xlab("Lake Turkana cichlids") + ylab("Sudanese cichlids") +
  theme(legend.position = "none",
        plot.tag = element_text(face = "bold"),
        plot.tag.position = c(0.01,1)) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  geom_abline(
    mapping = NULL,
    data = NULL,
    intercept = 0, slope=1,
    na.rm = FALSE,
    show.legend = NA, colour="white"
  )+ labs(tag = "I") 

grid.arrange(g1, g2, g3, g4, g5, g6, ncol=3)



#### Fig. S8: Victoria vs Southern Generalists for Western Lakes introgression #####
library(Hmisc)
dstats<-read.table("D:/Dropbox/victoriaGenomes/Dstats/Dstats.nubila.results",header=T)

# all numbers are negative, switch P1 and P2
dstats$extra<-dstats$P1
dstats$P1<-dstats$P2
dstats$P2<-dstats$extra
dstats$D<-(-dstats$D)
dstats$z<-(-dstats$z)
dstats$extra<-dstats$ABBA
dstats$ABBA<-dstats$BABA
dstats$BABA<-dstats$extra

dstats$low<-dstats$D-3*dstats$sd
dstats$high<-dstats$D+3*dstats$sd

nubsEd<-dstats[dstats$P3=="Edward",]
nubsKi<-dstats[dstats$P3=="Kivu",]
nubsAl<-dstats[dstats$P3=="Albert",]
nubsSa<-dstats[dstats$P3=="Saka",]
nubsTu<-dstats[dstats$P3=="Turkana",]
nubsSu<-dstats[dstats$P3=="Sudan",]
nubsEg<-dstats[dstats$P3=="Egypt",]

par(mfrow=c(1,1),mar=c(2,4,1,1),xaxs="r",yaxs="i")

plot(c(nubsEd$D[1],NA,nubsEd$D[-1],NA,NA,NA,NA,nubsKi$D[1],NA,nubsKi$D[-1],
       NA,NA,NA,NA,nubsAl$D[1],NA,nubsAl$D[-1],NA,NA,NA,NA,nubsSa$D[1],NA,nubsSa$D[-1],
       NA,NA,NA,NA,nubsEg$D[1],NA,nubsEg$D[-1],NA,NA,NA,NA,nubsSu$D[1],NA,nubsSu$D[-1],
       NA,NA,NA,NA,nubsTu$D[1],NA,nubsTu$D[-1]),
     ylim=c(0,0.1),cex=0.5,ylab="D",las=2,xlab="",xaxt="n")
errbar(x=1:length(c(nubsEd$D[1],NA,nubsEd$D[-1],NA,NA,NA,NA,nubsKi$D[1],NA,
                    nubsKi$D[-1],NA,NA,NA,NA,nubsAl$D[1],NA,nubsAl$D[-1],NA,NA,NA,NA,
                    nubsSa$D[1],NA,nubsSa$D[-1],NA,NA,NA,NA,nubsEg$D[1],NA,nubsEg$D[-1],
                    NA,NA,NA,NA,nubsSu$D[1],NA,nubsSu$D[-1],
                    NA,NA,NA,NA,nubsTu$D[1],NA,nubsTu$D[-1])),
       add=T,y=c(nubsEd$D[1],NA,nubsEd$D[-1],NA,NA,NA,NA,nubsKi$D[1],NA,nubsKi$D[-1],NA,NA,NA,NA,nubsAl$D[1],
                 NA,nubsAl$D[-1],NA,NA,NA,NA,nubsSa$D[1],NA,nubsSa$D[-1],NA,NA,NA,NA,nubsEg$D[1],NA,nubsEg$D[-1],
                 NA,NA,NA,NA,nubsSu$D[1],NA,nubsSu$D[-1],NA,NA,NA,NA,nubsTu$D[1],NA,nubsTu$D[-1]),
       yminus=c(nubsEd$low[1],NA,nubsEd$low[-1],NA,NA,NA,NA,nubsKi$low[1],NA,nubsKi$low[-1],NA,NA,NA,NA,
                nubsAl$low[1],NA,nubsAl$low[-1],NA,NA,NA,NA,nubsSa$low[1],NA,nubsSa$low[-1],
                NA,NA,NA,NA,nubsEg$low[1],NA,nubsEg$low[-1],NA,NA,NA,NA,nubsSu$low[1],NA,
                nubsSu$low[-1],NA,NA,NA,NA,nubsTu$low[1],NA,nubsTu$low[-1]),
       yplus=c(nubsEd$high[1],NA,nubsEd$high[-1],NA,NA,NA,NA,nubsKi$high[1],NA,nubsKi$high[-1],
               NA,NA,NA,NA,nubsAl$high[1],NA,nubsAl$high[-1],NA,NA,NA,NA,nubsSa$high[1],NA,nubsSa$high[-1],
               NA,NA,NA,NA,nubsEg$high[1],NA,nubsEg$high[-1],NA,NA,NA,NA,nubsSu$high[1],NA,nubsSu$high[-1],
               NA,NA,NA,NA,nubsTu$high[1],NA,nubsTu$high[-1]),
cap=0,col=NA,pch="",lwd=2.2,errbar.col=c(rep("blue",times=44),rep("darkblue",times=44),
                                          rep("cornflowerblue",times=44),rep("deepskyblue",times=44),
                                          rep("#1DA180",times=44),rep("#1DA180",times=44),rep("#1DA180",times=44)))
points(c(nubsEd$D[1],NA,nubsEd$D[-1],NA,NA,NA,NA,nubsKi$D[1],NA,nubsKi$D[-1],NA,NA,NA,NA,nubsAl$D[1],NA,
         nubsAl$D[-1],NA,NA,NA,NA,nubsSa$D[1],NA,nubsSa$D[-1],NA,NA,NA,NA,nubsEg$D[1],NA,nubsEg$D[-1],
         NA,NA,NA,NA,nubsSu$D[1],NA,nubsSu$D[-1],NA,NA,NA,NA,nubsTu$D[1],NA,nubsTu$D[-1]),
       cex=rep(c(1,1,1,rep(0.5,times=41)),times=7),
       pch=rep(c(17,rep(19,times=43)),times=7))
axis(1,at=c(22,43+22,43+43+22,43+43+43+22,43+43+43+43+22,43+43+43+43+43+22,43+43+43+43+43+43+22),
     labels = c("Edward","Kivu","Albert","Saka","Egypt","Sudan","Turkana"))
legend("topleft",pch=c(17,19,20),cex=c(1,1,1),title="P3:",title.adj = 0.05,
       legend = c("Northern Generalists","A. sp. \"nubila rocks\"","other Victoria Radiation"),bty="n")


