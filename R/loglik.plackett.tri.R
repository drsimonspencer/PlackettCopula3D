#' Log-likelihood for a trivariate Plackett copula
#'
#' Computes the log-likelihood of observed counts in a 3-dimensional binary
#' table under a trivariate Plackett copula model. Returns `-Inf` if any
#' computed probability is negative or `NA`.
#'
#' @param params named vector of parameters of the trivariate Plackett copula
#'  Must contain:
#'   \itemize{
#'     \item \code{u, v, w} : marginal prevalences for the three binary variables
#'     \item \code{phi.UV, phi.UW, phi.VW} : bivariate Plackett copula parameters
#'       for the corresponding pairs
#'     \item \code{psi} : dependence parameter for trivariate copula
#'   }
#' @param counts numeric vector of length 8 with observed counts for all
#'   combinations of the three binary variables (n000, n001, ..., n111)
#'
#' @returns numeric, the log-likelihood value
#' @export
#'
#' @examples
#' params <- c(u = 0.3, v = 0.4, w = 0.2, phi.UV = 2, phi.UW = 1.5, phi.VW = 3, psi=1.2)
#' counts <- c(n000=5, n001=3, n010=2, n011=1, n100=4, n101=2, n110=3, n111=6)
#' loglik.plackett.tri(params = params, counts = counts)

loglik.plackett.tri <- function(params, counts) {
  # Get counts in the same order as probs
  counts <- counts[c("n000","n001","n010","n011","n100","n101","n110","n111")]
  if (length(counts)!=8 | any(is.na(counts))) {stop("counts must have correct labels.\n")}

  probs <- prob.fun.tri(params)

  wh<-which(counts>0)
  if (any(is.na(probs[wh])) || any(probs[wh] < 0)) {
    #cat("Likelihood error:",probs,";",params,"\n")
    return(-Inf)
  } else {
    return(sum(counts[wh] * log(probs[wh])))
  }
}
