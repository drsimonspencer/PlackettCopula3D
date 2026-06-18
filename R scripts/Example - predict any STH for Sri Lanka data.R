# Script to compare anySTH predictions against the truth for Sri Lanka data
library(PlackettCopula3D)
data("SriLankaSTH")

# Calculate sample size and prevalences
samples<-apply(SriLankaSTH[,c("n000","n001","n010","n011","n100","n101","n110","n111")],1,sum)
ascaris<-apply(SriLankaSTH[,c("n100","n101","n110","n111")],1,sum)/samples
hookworm<-apply(SriLankaSTH[,c("n010","n011","n110","n111")],1,sum)/samples
trichuris<-apply(SriLankaSTH[,c("n001","n011","n110","n111")],1,sum)/samples

# Make predictions and print them
predictions<-predict.helminth(ascaris,hookworm,trichuris)
predictions$anySTH

# Plot predictions against the observations
plot(0:1,0:1,t="l",col="grey",xlim=c(0,1),ylim=c(0,1),xlab="Observed anySTH",ylab="predicted anySTH",main="Sri Lanka data")
points(1-SriLankaSTH[,c("n000")]/samples,predictions$anySTH)
