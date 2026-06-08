#' Sample from a multinomial distribution, conditional on its margins.
#'
#' Draw a random sample from a multinomial with probabilities given by the trivariate Plackett
#' copula, conditional on the known margins.
#'
#' @details
#' Samples are drawn in sequence, and dimensions are deleted if a draw from that dimension
#' would violate the marginal constraints.
#' @param size Number of samples to draw
#' @param probs Vector of length 8 giving the probabilities
#' @param margins Integer vector of length 3 giving the constraints.
#'
#' @returns Vector of length 8, giving the number of samples in each state.
#'
rmultinomPlackett<-function(size,probs,margins) {
  n.0<-size-margins
  n.1<-margins
  out<-rep(0,8)
  p<-probs
  p[which(p<0)]<-0
  it<-0
  while (max(n.0,n.1)>0) {
    it<-it+1
    # Remove categories that violate constraints
    if (n.0[1]==0) {p[c("p000","p001","p010","p011")]<-0}
    if (n.1[1]==0) {p[c("p100","p101","p110","p111")]<-0}
    if (n.0[2]==0) {p[c("p000","p001","p100","p101")]<-0}
    if (n.1[2]==0) {p[c("p010","p011","p110","p111")]<-0}
    if (n.0[3]==0) {p[c("p000","p010","p100","p110")]<-0}
    if (n.1[3]==0) {p[c("p001","p011","p101","p111")]<-0}
    # renormalise
    if (sum(p)>0) {
      p<-p/sum(p)
    } else {
      cat("Post pred sampling failure (retrying):",n.0,n.1,";",probs,"\n")
      return(rmultinomPlackett(size,probs,margins))
    }
    # Calculate how many to hit next margin
    n_sample<-min(n.0[which(n.0>0)],n.1[which(n.1>0)])
    s<-drop(rmultinom(1,n_sample,p))
    n.0[1]<-n.0[1]-sum(s[c("p000","p001","p010","p011")])
    n.1[1]<-n.1[1]-sum(s[c("p100","p101","p110","p111")])
    n.0[2]<-n.0[2]-sum(s[c("p000","p001","p100","p101")])
    n.1[2]<-n.1[2]-sum(s[c("p010","p011","p110","p111")])
    n.0[3]<-n.0[3]-sum(s[c("p000","p010","p100","p110")])
    n.1[3]<-n.1[3]-sum(s[c("p001","p011","p101","p111")])
    out<-out+s
  }
  #cat(sum(out),"conditional samples in",it,"iterations.\n")
  names(out)<-c("n000","n001","n010","n011","n100","n101","n110","n111")
  return(out)
}
