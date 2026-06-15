location<-"full" # "full", "randomized"
pdf_output<-FALSE
load(paste0("output/",location,"MCMC.Rdata"))

data("randomizedSTH") # For meta data
M<-nrow(randomizedSTH)
cols<-rep("red",length(randomizedSTH$CONTINENT))
cols[which(randomizedSTH$CONTINENT=="Asia")]<-"blue"
cols[which(randomizedSTH$CONTINENT=="South America")]<-"darkgreen"
pchs<-rep(1,length(randomizedSTH$CONTINENT))
pchs[which(randomizedSTH$BASELINE==0)]<-4
pchs[c(36,40,44,48,53,57,61,79,109,111,113,127)]<-0 # Can't think of a way to make this more general!

# Extract a subset? In the paper we show BASELINE and
STUDY<-substr(rownames(randomizedSTH),1,6)
study<-"MDA"# "SZZJSA","SYXMFQ","SWPPBC","STYOJP","STSVZC","STASBQ","SIZIIK","SHSFBR","SGNKUY","SEUXMF","SCNWRH","SALDNZ","SKA","BASELINE","MDA"
studylabel<-study
plotlegend<-TRUE
if (study=="") {
  wh<-1:length(randomizedSTH$BASELINE) # No subset
} else if (study=="BASELINE") {
  wh<-which(pchs==1)
  studylabel<-"Baseline surveys"
  plotlegend<-FALSE
} else if (study=="MDA") {
  wh<-which(pchs==4)
  studylabel<-"Post-treatment surveys"
  plotlegend<-FALSE
} else if (study=="PLACEBO") {
  wh<-which(pchs==0)
  studylabel<-"Post-placebo surveys"
  plotlegend<-FALSE
} else if (study=="SKA") {
  wh<-which(randomizedSTH$COUNTRY==study)
} else if (study=="PHL") {
  wh<-which(randomizedSTH$COUNTRY==study)
} else {
  wh<-which(STUDY==study)
}

# Make funnel plots
if (pdf_output) {
  pdf(paste0("output/BayesianFunnelPlots",location,study,".pdf"),width=6,height=6,pointsize=8)
} else {
  postscript(paste0("output/BayesianFunnelPlots",location,study,".eps"),width=6,height=6,pointsize=8)
}

par(mfrow=c(2,2))
m<-max(abs(intervals))
mains<-c("Ascaris-Hookworm interaction","Ascaris-Trichuris interaction","Hookworm-Trichuris interaction","Ascaris-Hookworm-Trichuris interaction")
labs<-c(expression(log(phi[UV])),expression(log(phi[UW])),expression(log(phi[VW])),expression(log(psi)))
for (j in 1:4) {
  set.seed(100)
  priors<-log(qf(c(0.025,0.5,0.975),mcmc$priors$d[j],mcmc$priors$d[j]))
  priors.sd<-sd(log(rf(100000,mcmc$priors$d[j],mcmc$priors$d[j])))
  s1<-seq(0,priors.sd,length.out=101)[2:101]
  m1<-s1^2*(s1^-2-priors.sd^-2)*meta.interval[2,j]
  lwr<-m1-qnorm(0.975)*s1
  upr<-m1+qnorm(0.975)*s1
  plot(0:1,t="n",main=mains[j],ylim=c(0,max(sds)),ylab="posterior sd",xlab=labs[j],xlim=c(-m,m))
  polygon(c(rev(lwr),meta.interval[2,j],upr),c(rev(s1),0,s1),
          density=NA,col="#DDDDDD",border=NA)
  abline(v=0,col="#BBBBBB")
  segments(intervals[wh,1,j],sds[wh,j],intervals[wh,3,j],sds[wh,j],col=cols[wh])
  points(intervals[wh,2,j],sds[wh,j],col=cols[wh],pch=pchs[wh])
  segments(meta.interval[1,j],meta.sd[j],meta.interval[3,j],meta.sd[j])
  points(meta.interval[2,j],meta.sd[j],pch=5)
  segments(priors[1],priors.sd,priors[3],priors.sd,col="darkorange")
  points(priors[2],priors.sd,col="darkorange",pch=5)
  legend("bottomleft",c("Prior","Asia","Africa","Brazil","All combined"),lty=1,col=c("darkorange","blue","red","darkgreen","black"),pch=c(5,1,1,1,5))
  if (plotlegend) {legend("bottomright",c("Baseline","Post-treatment","placebo"),pch=c(1,4,0))}
}
mtext(studylabel,line=-1.5,outer=TRUE)
dev.off()

# Make prevalence plots
if (pdf_output) {
  pdf(paste0("output/Prevalences",location,study,".pdf"),width=6,height=5,pointsize=8)
} else {
  postscript(paste0("output/Prevalences",location,study,".eps"),width=6,height=5,pointsize=8)
}
par(mfrow=c(1,3),mar=c(5, 2, 4, 2) + 0.1)
coco<-paste(randomizedSTH$CONTINENT,randomizedSTH$COUNTRY)
o<-order(coco[wh])
mains<-paste(c("Ascaris","Hookworm","Trichuris"),"prevalence")
xlabs<-c("u","v","w")
ylabs<-randomizedSTH$COUNTRY[wh[o]]
yax<-length(wh):1
blank<-c()
for (i in length(ylabs):2) {
  if (ylabs[i]==ylabs[i-1]) {
    ylabs[i]<-NA
  } else {
    yax[1:(i-1)]<-yax[1:(i-1)]+1
    blank<-c(blank,yax[i]+1)
  }
}
for (j in 1:3) {
  qnts<-apply(mcmc$chains[keep,(j-1)*M+wh],2,quantile,probs=c(0.025,0.5,0.975))
  plot(qnts[2,o],yax,col=cols[wh[o]],pch=pchs[wh[o]],yaxt="n",xlab=xlabs[j],xlim=c(0,1),main=mains[j],ylab="",ylim=c(4,max(yax)-3))
  axis(2,yax,ylabs,las=1,tick=FALSE,line=-0.75)
  segments(0,blank,1,blank,col="#DDDDDD",lty=3)
  segments(qnts[1,o],yax,qnts[3,o],yax,col=cols[wh[o]])
  if (j==1) {legend("topright",c("baseline","post-treatment","placebo"),pch=c(1,4,0),bg="white")}
}
dev.off()

# Plot posterior histograms
if (pdf_output) {
  pdf(paste0("output/Posteriors",location,".pdf"),width=6,height=6,pointsize=8)
} else {
  postscript(paste0("output/Posteriors",location,".eps"),width=6,height=6,pointsize=8)
}
par(mfrow=c(2,2),mar=c(5, 4, 4, 2) + 0.1)
plot.posteriors(mcmc,burnin=keep[1]-1,legs=c(rep("topleft",3),"topright"),ylim=c(0,3.6))
dev.off()

