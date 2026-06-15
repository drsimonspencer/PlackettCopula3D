library(PlackettCopula3D)

# load the randomized data set.
if (!exists("d")) {
  data("randomizedSTH")
  d<-randomizedSTH
  randomized<-TRUE
} else {
  randomized<-FALSE
}

# Perform MCMC analysis of the complete data
M<-nrow(d)
burnin<-100
N<-10000+burnin
keep<-(burnin+1):N
start_time<-Sys.time()
mcmc<-mcmc.plackett(d,N=N) # (takes about 30 minutes on my laptop)
end_time<-Sys.time()
print(difftime(end_time,start_time,units="mins"))

# Plot MCMC diagnostics
plot.mcmc.diagnostics(mcmc,"trace")

# Save posterior
if (randomized) {
  save(mcmc,file="output/posterior_randomized.rda")
} else {
  mcmc$counts[,]<-NA
  save(mcmc,file="data/posterior.rda")
}

# Calculate 95% Credible intervals and posterior sd
meta.interval<-apply(log(mcmc$chains[keep,(3*M+1):(3*M+4)]),2,quantile,probs=c(0.025,0.5,0.975))
meta.sd<-apply(log(mcmc$chains[keep,(3*M+1):(3*M+4)]),2,sd)

# Compare phi and psi by analysing each dataset individually
mcmc.list<-list()
intervals<-array(NA,c(M,3,4))
sds<-matrix(NA,M,4)
dimnames(intervals)<-list(rownames(d),c("lwr","med","upr"),c("phi_UV","phi_UW","phi_VW","psi"))
for (i in 1:M) { # Analyse each dataset individually (takes some time)
  cat("Survey",i,"\n")
  mcmc.list[[i]]<-mcmc.plackett(d[i,],N=N)
  intervals[i,,]<-apply(log(mcmc.list[[i]]$chains[keep,4:7]),2,quantile,probs=c(0.025,0.5,0.975))
  sds[i,]<-apply(log(mcmc.list[[i]]$chains[keep,4:7]),2,sd)
  par(mfrow=c(2,4))
  plot.mcmc.diagnostics(mcmc.list[[i]],"trace")
  mtext(paste("Survey",i),outer=TRUE,line=-1)
}

Sample.size <- rowSums(d[,c("n000","n001","n010","n011","n100","n101", "n110", "n111")])
truth<-1-d$n000/Sample.size

# Store key output in output folder
if (randomized) {
  save(mcmc,keep,truth,meta.interval,meta.sd,intervals,sds,file="output/randomizedMCMC.Rdata")
} else {
  save(mcmc,keep,truth,meta.interval,meta.sd,intervals,sds,file="output/fullMCMC.Rdata")
}

