location<-"full" # "full", "randomized"
load(paste0("output/",location,"MCMC.Rdata"))
library(PlackettCopula3D)

# load International data from classic paper
data("InternationalSTH")
InternationalSTH$anySTHindep<-1-(1-InternationalSTH$Ascaris)*(1-InternationalSTH$Hookworm)*(1-InternationalSTH$Trichuris)

# Deal with our data
data("randomizedSTH") # for meta data and marginal prevalences
Sample.size <- rowSums(randomizedSTH[,c("n000","n001","n010","n011","n100","n101", "n110", "n111")])
M<-length(randomizedSTH$COUNTRY)
sum.1 <- cbind(rowSums(randomizedSTH[,c("n100", "n101", "n110", "n111")]),
               rowSums(randomizedSTH[,c("n010", "n011", "n110", "n111")]),
               rowSums(randomizedSTH[,c("n001", "n101", "n011", "n111")]))
colnames(sum.1)<-c("Ascaris","Hookworm","Trichuris")

# Fit and predict for classic model
lm0<-lm(Any.worm.observed~0+anySTHindep,data=InternationalSTH)
predictions<-predict(lm0,list(anySTHindep=1-apply(1-sum.1/Sample.size,1,prod)))

# Make predictions for our model for our data
pred.list<-list()
ppintervals<-matrix(NA,M,3)
dimnames(ppintervals)<-list(rownames(randomizedSTH),c("lwr","med","upr"))
mxintervals<-matrix(NA,M,3)
set.seed(101)
for (i in 1:M) { # Predict for each dataset
  cat("Location",i,"\n")
  mxintervals[i,]<-c(max(sum.1[i,1],sum.1[i,2],sum.1[i,3])/Sample.size[i],NA,min(1,sum(sum.1[i,])/Sample.size[i]))
  # Now predict the any STH data
  pred.list[[i]]<-predict.plackett(mcmc,prevalences=sum.1[i,]/Sample.size[i],sizes=rep(Sample.size[i],3))
  ppintervals[i,]<-(Sample.size[i]-quantile(pred.list[[i]]$draws[keep,1],probs=c(0.025,0.5,0.975)))/Sample.size[i]
}

RMSE<-sqrt(c(mean((ppintervals[,2]-truth)^2),mean((predictions-truth)^2)))
RMSRE<-sqrt(c(mean(((ppintervals[,2]-truth)/truth)^2),mean(((predictions-truth)/truth)^2)))
copulaIsWorse<-(ppintervals[,2]-truth)^2>(predictions-truth)^2
copulaIsBetter<-(ppintervals[,2]-truth)^2<(predictions-truth)^2
outOfBounds<-(predictions<mxintervals[,1]) | (predictions>mxintervals[,3])

save(ppintervals,mxintervals,truth,RMSE,RMSRE,outOfBounds,file=paste0("output/",location,"predictions.Rdata"))


