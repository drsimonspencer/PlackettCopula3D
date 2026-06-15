location<-"full" # "full", "randomized"
load(paste0("output/",location,"MCMC.Rdata"))
library(PlackettCopula3D)

# load International data from classic paper
data("InternationalSTH")
attach(InternationalSTH)
M<-length(Sample.size)
InternationalSTH$anySTHindep<-1-(1-Ascaris)*(1-Hookworm)*(1-Trichuris)

# Fit and predict for classic model
lm0<-lm(Any.worm.observed~0+anySTHindep,data=InternationalSTH)
predictions<-predict(lm0)

# Make predictions for our model for classic data
pred.list<-list()
ppintervals<-matrix(NA,M,3)
dimnames(ppintervals)<-list(Community,c("lwr","med","upr"))
mxintervals<-matrix(NA,M,3)
set.seed(101)
for (i in 1:M) { # Analyse each dataset individually
  cat("Location",i,"\n")
  mxintervals[i,]<-c(max(Ascaris[i],Hookworm[i],Trichuris[i]),NA,min(1,Ascaris[i]+Hookworm[i]+Trichuris[i]))
  # Now predict the any STH data
  pred.list[[i]]<-predict.plackett(mcmc,prevalences=c(Ascaris[i],Hookworm[i],Trichuris[i]),sizes=rep(Sample.size[i],3))
  ppintervals[i,]<-(Sample.size[i]-quantile(pred.list[[i]]$draws[keep,1],probs=c(0.025,0.5,0.975)))/Sample.size[i]
}

RMSE<-sqrt(c(mean((ppintervals[,2]-Any.worm.observed)^2),mean((predictions-Any.worm.observed)^2)))
RMSRE<-sqrt(c(mean(((ppintervals[,2]-Any.worm.observed)/Any.worm.observed)^2),mean(((predictions-Any.worm.observed)/Any.worm.observed)^2)))
copulaIsWorse<-(ppintervals[,2]-Any.worm.observed)^2>(predictions-Any.worm.observed)^2
copulaIsBetter<-(ppintervals[,2]-Any.worm.observed)^2<(predictions-Any.worm.observed)^2
outOfBounds<-(predictions<mxintervals[,1]) | (predictions>mxintervals[,3])

save(ppintervals,mxintervals,Any.worm.observed,RMSE,RMSRE,outOfBounds,file=paste0("output/",location,"predictionsInternational.Rdata"))
detach(InternationalSTH)
