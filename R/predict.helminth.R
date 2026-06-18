#' Make predictions for new locations for ascaris, hookworm and trichuris
#'
#' This function takes MCMC output and extracts the Plackett copula parameters. It then predicts
#' the prevalence of each of the 8 states based on marginal prevalence data from a new location.
#' Data in the new locations can be of the form of survey counts, or as a fixed (known) prevalence.
#'
#' @param ascaris Ascaris prevalence vector (with entries between 0 and 1).
#' @param hookworm Hookworm prevalence vector (with entries between 0 and 1).
#' @param trichuris Trichuris prevalence vector (with entries between 0 and 1).
#' @param confidence confidence level for the credible interval, default 0.95.
#'
#' All three vectors must ha
#'
#' @returns A list with:
#' \item{anySTH}{Predictions for the prevalence of each of the 8 possible states, plus "any STH".}
#' \item{credibleInterval}{Approximate credible interval for anySTH.}
#' \item{predictions}{Predictions for the prevalence of each of the 8 possible states, plus "any STH".}
#'
#' @examples
#' predict.helminth(ascaris=0.5,hookworm=0.6,trichuris=0.2)
#' predict.helminth(ascaris=c(0.1,0.2),hookworm=c(0.1,0.3),trichuris=c(0.2,0.1))
#'
#' @export
predict.helminth<-function(ascaris, hookworm, trichuris,confidence=0.95) {
  if (length(ascaris)!=length(hookworm) || length(hookworm)!=length(trichuris)) {stop("ascaris, hookworm and trichuris must have the same length.\n")}
  if (max(ascaris,hookworm,trichuris)>1 | min(ascaris,hookworm,trichuris)<0) {stop("prevalences must be between 0 and 1.\n")}
  data("thinned_posterior")
  quantiles<-apply(thinned_posterior$chains[,c("phi.UV","phi.UW","phi.VW","psi")],2,quantile,probs=c((1-confidence)/2,0.5,0.5+confidence/2))
  N<-length(ascaris)
  predictions<-matrix(NA,nrow = N, ncol=8)
  colnames(predictions)<-c("p000","p001","p010","p011","p100","p101","p110","p111")
  credibleInterval<-matrix(NA,nrow = N, ncol=2)
  colnames(credibleInterval)<-c("lwr","upr")
  anySTH<-rep(NA,N)
  for (i in 1:N) {
    #phi.UV   phi.UW   phi.VW       psi
    #2.5%  1.922374 2.635057 2.125477 0.7309506
    #50%   2.117274 2.905510 2.322982 0.9317213
    #97.5% 2.340432 3.212645 2.543451 1.1894187
    params<-c(ascaris[i],hookworm[i],trichuris[i],quantiles[2,])
    names(params)[1:3]<-c("u","v","w")
    params_lwr<-c(params[1:3],quantiles[3,1:3],psi=quantiles[1,4])
    params_upr<-c(params[1:3],quantiles[1,1:3],psi=quantiles[3,4])
    predictions[i,]<-prob.fun.tri(params)
    anySTH[i]<-1-predictions[i,1]
    credibleInterval[i,]<-c(1-prob.fun.tri(params_lwr)[1],1-prob.fun.tri(params_upr)[1])
  }
  names(anySTH)<-NULL
  return(list(anySTH=anySTH,credibleInterval=credibleInterval,predictions=predictions))
}
