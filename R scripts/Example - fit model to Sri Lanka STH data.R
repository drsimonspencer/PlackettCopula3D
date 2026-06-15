# Script to analyse the Sri Lanka STH dataset

library(PlackettCopula3D)
data("SriLankaSTH")

# Check dataset has columns called n000, n001, n010, n011, n100, n101, n110 and n111.
# The order of the digits is Ascaris, Hookworm and Trichuris in this dataset.
SriLankaSTH

burnin<-100
N<-2000+burnin
# Run MCMC
mcmc<-mcmc.plackett(SriLankaSTH,N=N)

# Plot MCMC diagnostics
plot.mcmc.diagnostics(mcmc,"trace")

# Plot posteriors
par(mfrow=c(2,2))
plot.posteriors(mcmc,burnin=burnin,legs=c(rep("topleft",3),"topright"),ylim=c(0,3.6))

