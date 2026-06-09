#' Make predictions for new locations from Plackett Copula Parameters
#'
#' This function takes MCMC output and extracts the Plackett copula parameters. It then predicts
#' the prevalence of each of the 8 states based on marginal prevalence data from a new location.
#' Data in the new locations can be of the form of survey counts, or as a fixed (known) prevalence.
#'
#' @param mcmc Output list from mcmc.plackett().
#' @param prevalences Vector of length 3, giving the prevalences u, v and w.
#' @param sizes Vector of length 3 giving the sample size (number of people surveyed).
#' If prevalences are to be treated as fixed, then leave as NULL.
#' @param priors List containing the parameters of the Beta prior for \eqn{u}, \eqn{v}, and \eqn{w}.
#'   Entries must include a.u, b.u, a.v, b.v, a.w, b.w (default: 0.5). Not used if size is NULL.
#' @param max.attempts Maximum number of attempts to simulate (default 100).
#'
#' @returns A list with:
#' \item{predictions}{Predictions for the prevalence of each of the 8 possible states, plus "any STH".}
#' \item{params}{Parameters for the prevalence of each of the 8 possible states, plus "any STH".}
#' \item{draws}{Samples from the posterior predictive distribution (requires sizes to be equal).}
#' \item{prevalences}{Vector of prevalences.}
#' \item{sizes}{Vector of sizes.}
#' \item{priors}{List of priors used.}
#'
#' @examples
#' counts <- c(10, 5, 8, 6, 7, 4, 9, 11)
#' names(counts) <- c("n000","n001","n010","n011","n100","n101","n110","n111")
#' mcmc <- mcmc.plackett(counts, N = 1000)
#' predict.plackett(mcmc, prevalences=c(u=5/10,v=6/12,w=1/14), sizes=c(u=10,v=12,w=14))
#'
#' @export
predict.plackett<-function(mcmc, prevalences, sizes=NULL, priors=NULL, max.attempts=100) {
  if (length(prevalences)!=3) {stop("prevalences must be of length 3.\n")}
  if (max(prevalences)>1 | min(prevalences)<0) {stop("prevalences must be between 0 and 1.\n")}
  if (!is.null(sizes) && length(sizes)!=3) {stop("sizes must be of length 3.\n")}
  # MCMC iterations
  N<-dim(mcmc$chains)[1]
  # Number of previous locations
  M<-(dim(mcmc$chains)[2]-4)/3
  # Set default priors if not specified.
  if (is.null(priors$a.u)) {priors$a.u <- 0.5}
  if (is.null(priors$b.u)) {priors$b.u <- 0.5}
  if (is.null(priors$a.v)) {priors$a.v <- 0.5}
  if (is.null(priors$b.v)) {priors$b.v <- 0.5}
  if (is.null(priors$a.w)) {priors$a.w <- 0.5}
  if (is.null(priors$b.w)) {priors$b.w <- 0.5}
  predictions<-matrix(NA,nrow = N, ncol=9)
  colnames(predictions)<-c("p000","p001","p010","p011","p100","p101","p110","p111","p_any")
  params<-matrix(NA, nrow = N, ncol=7)
  colnames(params)<-c("u","v","w",colnames(mcmc$chains)[(3*M+1):(3*M+4)])
  draws<-matrix(NA,nrow = N, ncol=8)
  colnames(draws)<-c("n000","n001","n010","n011","n100","n101","n110","n111")
  # Loop through MCMC iterations
  for (i in 1:N) {
    attempt <- 0
    while (attempt<max.attempts) {
      attempt <- attempt + 1
      if (!is.null(sizes)) {
        a<-c(priors$a.u+sizes[1]*prevalences[1],
             priors$a.v+sizes[2]*prevalences[2],
             priors$a.w+sizes[3]*prevalences[3])
        b<-c(priors$b.u+sizes[1]*(1-prevalences[1]),
             priors$b.v+sizes[2]*(1-prevalences[2]),
             priors$b.w+sizes[3]*(1-prevalences[3]))
        p<-rbeta(3,a,b)
      } else {
        attempt<-max.attempts
        p<-prevalences
      }
      names(p)<-c("u","v","w")
      probs <- prob.fun.tri(c(p,mcmc$chains[i,(3*M+1):(3*M+4)]))
      if (!any(is.na(probs)) || any(probs<0)) {
        if (attempt>1) {cat("Attempts:",attempt,"\n")}
        attempt<-max.attempts
      }
    }
    predictions[i,]<-c(probs,1-probs["p000"])
    params[i,]<-c(p,mcmc$chains[i,(3*M+1):(3*M+4)])
    if (sizes[1]==sizes[2] && sizes[2]==sizes[3]) {
      draws[i,]<-rmultinomPlackett(sizes[1],probs,round(sizes[1]*prevalences,0))
    }
  }
  if (sizes[1]!=sizes[2] || sizes[2]!=sizes[3]) {draws<-NULL}
  return(list(predictions=predictions,params=params,draws=draws,prevalences=prevalences,sizes=sizes,priors=priors))
}
