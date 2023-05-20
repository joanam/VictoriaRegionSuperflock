require(tidyr)
f4lv<-read.table("D:/Dropbox/victoriaGenomes/Dstats/LV.4pop.ecoGroups.results.noLakeOrder.txt",header=T)

# Function that reformats the f4 values for plotting
reformatLV<-function(f4){
  # Add lake information and sort by lake
  f4<-separate(data=f4,col=P1,into = "lake",sep = "_",extra = "drop",remove=F)
  f4<-separate(data=f4,col=P3,into = c("tmp","eco1"),sep = "_",extra = "merge",remove=F)
  f4<-separate(data=f4,col=P4,into = c("tmp","eco2"),sep = "_",extra = "merge",remove=F)
  f4$lake<-factor(f4$lake,c("Edward","Kivu","Albert"))
  f4$col<-f4$lake
  levels(f4$col)<-c("blue","darkblue","cornflowerblue")
  
  # Add a blank line after each lake group
  f4$sd<-as.double(f4$sd)
  f4$f4<-as.double(f4$f4)
  f4$z<-as.double(f4$z)
  f4$lake<-factor(f4$lake,c("Edward","Kivu","Albert"))
  f4$col<-factor(f4$col,c("blue","darkblue","cornflowerblue"))
  return(f4)
}

f4lv<-reformatLV(f4lv)

#### Fig. 5B: Plot of f4 values   ####################################################################


# Plot f4 full dataset

par(mar=c(2,5,1,1),xaxs="i",cex.lab=1.5,mgp=c(2,0.5,0))
plot(f4lv$f4,col=NA,ylim=c(-0.9e-4,3.8e-4),xlab="",ylab="f4",
     xaxt="n",yaxt="n",xlim=c(0,length(f4lv$col)+1))

#  Plot error bars (black if significant)
arrows(x0=1:length(f4lv$P1), y0=f4lv$f4-f4lv$sd*3, 
       x1=1:length(f4lv$P1), y1=f4lv$f4+f4lv$sd*3, 
       code=3, angle=90, length=0.03, 
       col=ifelse(f4lv$z>3,"black","grey"), lwd=2)

# Add a horizontal line at 0
abline(h=0)

# Add vertical lines that separate the ecomorph pairs
abline(v=which(is.na(f4lv$total)),col="darkgrey")

# Add y axis
axis(2,at=c(0,-1.5e-4,-3e-4,3e-4,1.5e-4))

# Label the ecomorph pairs
text(x = c(1:length(f4lv$eco2)),
     y = -0.00007,labels = f4lv$Sym1,
     xpd=T,cex=1.5)
text(x = c(1:length(f4lv$eco2)),
     y = -0.000095,labels = f4lv$Sym2,
     xpd=T,cex=1.5)

# Add points in lake colour for the f4 values
points(f4lv$f4,col=as.character(f4lv$col),pch=19,cex=1.2)

# Add legends
legend("topleft",legend = c("Albert","Edward","Kivu"),
       horiz = F,box.lty=0,pch=19,title="Victoria versus",
       col=c("cornflowerblue","blue","darkblue"),cex=1.4,inset=0.01)
legend("topright",inset=0.01,
       legend=c("in   insectivore","pi   piscivore","el   epilithic algae scraper",
                "ep  epiphylic algae scraper","de  detritivore",
                "os  oral sheller snail eater","pc  pharyngeal snail crusher",
                "pa  paedopaghe"),ncol=2,cex=1.4,box.lty=0
)

