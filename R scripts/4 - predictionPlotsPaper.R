location<-"full" # "full", "randomized"
pdf_output<-TRUE
library(PlackettCopula3D)
#
# Open file
#
if (pdf_output) {
  pdf("output/PredictionPlot.pdf",width=6,height=3.5,pointsize=8)
} else {
  postscript("output/PredictionPlot.pdf",width=6,height=3.5,pointsize=8)
}
par(mfrow=c(1,2))
#
# Load international predictions
#
data("InternationalSTH")
load(paste0("output/",location,"predictionsInternational.Rdata"))
cols<-rep("red",length(InternationalSTH$Continent))
cols[which(InternationalSTH$Continent=="Asia")]<-"blue"
cols[which(InternationalSTH$Continent=="Latin America")]<-"darkgreen"
cols[which(InternationalSTH$Continent=="Oceania")]<-"brown"
#
# Plot posterior predictive against truth for international data
#
plot(c(0,1),c(0,1),t="l",col="lightgrey",xlim=c(0,1),ylim=c(0,1),
     xlab="observed proportion with any STH",ylab="prediction from marginal data",main="Data from previous study")
segments(Any.worm.observed,mxintervals[,1],Any.worm.observed,mxintervals[,3],lwd=2,col="lightblue")
segments(Any.worm.observed,ppintervals[,1],Any.worm.observed,ppintervals[,3],col=cols)
points(Any.worm.observed,ppintervals[,2],col=cols)
# points(truth,predictions,col="pink",pch=".")
#legend("bottomright","bounds",col="lightblue",lwd=2)
legend("topleft",c("Africa","Asia","Americas","Pacific"),col=c("red","blue","darkgreen","brown"),lty=1,pch=1)
#
# Load full data predictions
#
data("randomizedSTH")
load(paste0("output/",location,"predictions.Rdata"))
#
cols<-rep("red",length(randomizedSTH$CONTINENT))
cols[which(randomizedSTH$CONTINENT=="Asia")]<-"blue"
cols[which(randomizedSTH$CONTINENT=="South America")]<-"darkgreen"
pchs<-rep(1,length(randomizedSTH$CONTINENT))
pchs[which(randomizedSTH$BASELINE==0)]<-4
pchs[c(36,40,44,48,53,57,61,79,109,111,113,127)]<-0 # Can't think of a way to make this more general!
#
# Plot posterior predictive against truth for our data (location)
#
plot(c(0,1),c(0,1),t="l",col="lightgrey",xlim=c(0,1),ylim=c(0,1),
     xlab="observed proportion with any STH",ylab="prediction from marginal data",main="Data from this study")
segments(truth,mxintervals[,1],truth,mxintervals[,3],lwd=2,col="lightblue")
segments(truth,ppintervals[,1],truth,ppintervals[,3],col=cols)
points(truth,ppintervals[,2],col=cols,pch=pchs)
# points(truth,predictions,col="pink",pch=".")
legend("bottomright",c("Baseline","Post-treatment","Placebo"),pch=c(1,4,0))
legend("topleft","bounds",col="lightblue",lwd=2)
#
dev.off()

