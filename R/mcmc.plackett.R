#' MCMC for the Trivariate Plackett Copula
#'
#' Runs a Markov Chain Monte Carlo (MCMC) algorithm to estimate the parameters
#' of the trivariate Plackett copula from contingency tables of three binary variables.
#'
#' @details
#' The input `counts` must be a vector of length 8 containing the observed frequencies of the 8 possible
#' outcomes of three binary variables, with \code{names} as:
#' \code{n000, n001, n010, n011, n100, n101, n110, n111}. If more than one dataset is specified
#' `counts` must be a matrix with 8 columns and \code{colnames} as the same as above.
#'
#' Parameters
#' - \eqn{u, v, w}: marginal prevalences of each binary variable, dataset-specific.
#' - \eqn{\phi_{UV}, \phi_{UW}, \phi_{VW}}: bivariate dependence parameters.
#' - \eqn{\psi}: trivariate dependence parameter.
#'
#' Priors
#' - \eqn{u_s} ~ \eqn{\mathrm{Beta}(a_U^s, b_U^{(s)})}, for survey \eqn{s}, same for \eqn{v_s} and \eqn{w_s}.
#' Jeffrey's prior by default (\eqn{a_U^{(s)}=b_U^{(s)}=0.5}).
#' - \eqn{\phi_{UV}, \phi_{UW}, \phi_{VW}}, \eqn{\psi} ~ F(d, d), with \code{d = 3} by default.
#'
#' Proposals
#' - \eqn{u, v, w}: proposed from Beta distributions matching their marginal posterior's full conditionals,
#'   e.g. \eqn{\mathrm{Beta}(a^{U} + \text{counts for U=0}, b_{U} + \text{counts for U=1})}.
#'   Accepted/rejected with a Metropolis–Hastings step.
#'
#' - \eqn{\phi_{UV}, \phi_{UW}, \phi_{VW}}, \eqn{\psi}: proposed from log-Gaussian random walks,
#'   e.g. \eqn{log(\phi_{UV}')\sim\mathcal{N}(log(\phi_{UV}), \sigma^2)}.
#'   Proposal variances \eqn{\sigma^2} are adapted via Robbins–Monro updates to target
#'   an acceptance rate \code{mcmc.options$target.accept} (default 0.44).
#'
#' Multiple datasets
#' If `counts` has multiple rows (\eqn{M \times 8}), the algorithm updates \eqn{u, v, w}
#' separately for each dataset, while \eqn{\phi_{UV}, \phi_{UW}, \phi_{VW}}, and \eqn{\psi}
#' are shared across all datasets.
#'
#' If only one dataset is supplied, `counts` must be a named vector with \code{n000}, \code{n001}, \code{n010}, \code{n011}, \code{n100},
#'   \code{n101}, \code{n110}, \code{n111}.
#'
#' @param counts A matrix of dimension \eqn{M \times 8} (or vector of length 8) with counts
#'   of the binary outcomes.
#' @param priors List containing the parameters of the Beta prior for \eqn{u}, \eqn{v}, and \eqn{w}.
#'   Entries must include a.u, b.u, a.v, b.v, a.w, b.w (default: `rep(0.5,M)`).
#'   Degrees of freedom for the F-distributed priors for \eqn{\phi_{UV}}, \eqn{\phi_{UW}}, \eqn{\phi_{VW}} and \eqn{\psi} are in a vector
#'    of length 4 called `d` (default: `rep(3,4)`).
#' @param N Number of MCMC iterations (default: 5000).
#' @param mcmc.options A list containing
#'   `target.accept` - target acceptance rate for Robbins–Monro adaptation (default: 0.44).
#'   `delta` - Robbins–Monro adaptation constant for phi and psi (default: 1/(0.44*(1-0.44))).
#'   `delta.lambda` - Robbins-Monro adaptation contact for u, v and w (default: 0 = off).
#'
#' @returns A list with:
#' \item{chains}{Matrix of posterior samples for \eqn{u}, \eqn{v}, \eqn{w}, \eqn{\phi_{UV}},
#'   \eqn{\phi_{UW}}, \eqn{\phi_{VW}} and \eqn{\psi}.}
#' \item{p_any}{Matrix of posterior samples probability of any STH: \eqn{1-p_{000}}.}
#' \item{acceptance}{Acceptance rates for all parameters.}
#' \item{sigma.trace}{Trace of proposal SDs for \eqn{\phi_{UV}}, \eqn{\phi_{UW}}, \eqn{\phi_{VW}} and \eqn{\psi}.}
#' \item{lambda.trace}{Trace of adaptive parameters for \eqn{u}, \eqn{v}, \eqn{w}.}
#' \item{counts}{Input data of STH counts.}
#' \item{priors}{List of prior parameters.}
#'
#' @examples
#' # One dataset (must have colnames containing n000, n001, n010, n011, n100, n101, n110, n111)
#' counts <- c(10, 5, 8, 6, 7, 4, 9, 11)
#' names(counts) <- c("n000","n001","n010","n011","n100","n101","n110","n111")
#' out <- mcmc.plackett(counts, N = 1000)
#' cbind(out$chains.u,out$chains.v,out$chains.w,out$chains)
#'
#' # Analyse two datasets, but with shared interaction parameters.
#' counts.multi <- rbind(
#'   c(10, 5, 8, 6, 7, 4, 9, 11),
#'   c(15, 3, 6, 5, 10, 2, 8, 12)
#' )
#' colnames(counts.multi) <- names(counts)
#' out.multi <- mcmc.plackett(counts.multi, N = 1000)
#' out.multi$acceptance
#'
#' @export

mcmc.plackett <- function(counts,
                          priors = list(),
                          N = 5000,
                          mcmc.options = list(target.accept = 0.44, delta = 1/0.44/0.56, delta.lambda = 0)) {

  if (is.null(nrow(counts))) {
    counts <- matrix(counts, nrow = 1, dimnames = list("",names(counts)))
  } else if (is.data.frame(counts)) {
    counts <- as.matrix(counts[c("n000","n001","n010","n011","n100","n101","n110","n111")])
  }
  if (length(which(colnames(counts)%in%c("n000","n001","n010","n011","n100","n101","n110","n111")))!=8) {
    stop("Counts must have names or colnames n000,n001,n010,n011,n100,n101,n110,n111\n")
  }

  # Number of datasets
  M <- nrow(counts)

  if (is.null(rownames(counts))) {rownames(counts)<-1:M}

  # Set default priors if they are not specified
  if (is.null(priors$a.u)) {priors$a.u <- rep(0.5,M)}
  if (is.null(priors$b.u)) {priors$b.u <- rep(0.5,M)}
  if (is.null(priors$a.v)) {priors$a.v <- rep(0.5,M)}
  if (is.null(priors$b.v)) {priors$b.v <- rep(0.5,M)}
  if (is.null(priors$a.w)) {priors$a.w <- rep(0.5,M)}
  if (is.null(priors$b.w)) {priors$b.w <- rep(0.5,M)}
  if (is.null(priors$d)) {priors$d <- rep(3,4)}

  # Pre-compute sums for beta proposals
  sum.0 <- cbind(rowSums(counts[,c("n000", "n001", "n010", "n011"), drop = F]),
                 rowSums(counts[,c("n000", "n001", "n100", "n101"), drop = F]),
                 rowSums(counts[,c("n000", "n100", "n010", "n110"), drop = F]))
  sum.1 <- cbind(rowSums(counts[,c("n100", "n101", "n110", "n111"), drop = F]),
                 rowSums(counts[,c("n010", "n011", "n110", "n111"), drop = F]),
                 rowSums(counts[,c("n001", "n101", "n011", "n111"), drop = F]))

  # Initialize prevalences at Bayesian point estimates
  p <- cbind(c((priors$a.u+sum.1[,1])/(priors$a.u+priors$b.u+rowSums(counts)),use.names=FALSE),
             c((priors$a.v+sum.1[,2])/(priors$a.v+priors$b.v+rowSums(counts)),use.names=FALSE),
             c((priors$a.w+sum.1[,3])/(priors$a.w+priors$b.w+rowSums(counts)),use.names=FALSE))
  colnames(p) <- c("u","v","w")
  # Initialize copula parameters at 1 because MLEs maybe 0 or infinity
  params <- rep(1,4)
  names(params) <- c("phi.UV","phi.UW", "phi.VW", "psi")

  # Initialize sigmas and lambdas (MCMC proposal tuning parameters)
  sigma <- rep(1,4)
  names(sigma) <- names(params)
  lambda <- matrix(0, nrow = M, ncol = 3)
  dimnames(lambda) <- list(rownames(counts),c("u","v","w"))

  # Initialize storage
  chains <- matrix(NA, nrow = N, ncol = 3*M+4)
  colnames(chains) <- c(paste("u",rownames(counts)),paste("v",rownames(counts)),paste("w",rownames(counts)),names(params))
  p_any <- matrix(NA, nrow = N, ncol = M)
  colnames(p_any) <- rownames(counts)

  accept <- matrix(NA, nrow = N, ncol = 3*M+4)
  colnames(accept) <- colnames(chains)
  sigma.trace <- matrix(NA, nrow = N, ncol = 4)
  colnames(sigma.trace) <- names(params)
  lambda.trace <- matrix(NA, nrow = N, ncol = 3*M)
  colnames(lambda.trace) <- colnames(chains)[1:(3*M)]

  for(it in 1:N) {

    for(i in 1:M) {
      for (j in 1:3) {
        # ----- Update p
        proposal <- p[i,]
        proposal[j] <- rbeta(1,priors[[2*j-1]][i]+sum.1[i,j]+lambda[i,j]*p[i,j], priors[[2*j]][i]+sum.0[i,j]+lambda[i,j]*(1-p[i,j]))

        # log acceptance ratio: proposal ratio
        lar <- dbeta(p[i,j], priors[[2*j-1]][i]+sum.1[i,j]+lambda[i,j]*p[i,j], priors[[2*j]][i]+sum.0[i,j]+lambda[i,j]*(1-p[i,j]), log = TRUE) -
               dbeta(proposal[j],      priors[[2*j-1]][i]+sum.1[i,j]+lambda[i,j]*p[i,j], priors[[2*j]][i]+sum.0[i,j]+lambda[i,j]*(1-p[i,j]), log = TRUE)
        # prior ratio
        lar <- lar + dbeta(proposal[j], priors[[2*j-1]][i], priors[[2*j]][i], log = TRUE) -
                     dbeta(p[i,j]     , priors[[2*j-1]][i], priors[[2*j]][i], log = TRUE)
        # likelihood ratio
        lar <- lar + loglik.plackett.tri(c(proposal,params), counts[i,]) -
                     loglik.plackett.tri(c(p[i,]   ,params), counts[i,])

        # Accept wp exp(min(1,lar))
        if(!is.na(lar) && log(runif(1)) < lar){
          p[i,] <- proposal
          accept[it,M*(j-1)+i] <- 1
        } else {
          accept[it,M*(j-1)+i] <- 0
        }
        # --- Robbins-Monro update ---
        lambda[i,j] <- max(0,lambda[i,j] - (mcmc.options$delta.lambda) * (min(1,exp(lar)) - mcmc.options$target.accept))
        lambda.trace[it,M*(j-1)+i]<-lambda[i,j]
      }
    }
    # ------ Update phis and psi with RWM
    for (j in 1:4) {
      proposal <- params
      proposal[j] <- rlnorm(1,log(params[j]),sigma[j]) # log-normal proposal

      # Hastings ratio
      #lar <- dlnorm(params[j],log(proposal[j]),sigma[j],log=TRUE) - dlnorm(proposal[j],log(params[j]),sigma[j],log=TRUE)
      lar<-log(proposal[j])-log(params[j]) # Hastings ratio is just the Jacobians, it's symmetric on the log scale.

      # log acceptance ratio: priors
      lar <- lar + df(proposal[j], df1 = priors$d[j], df2 = priors$d[j], log = TRUE) - df(params[j], df1 = priors$d[j], df2 = priors$d[j], log = TRUE)


      # likelihood ratio
      for(i in 1:M) {
        lar <- lar + loglik.plackett.tri(params=c(p[i,],proposal), counts[i,]) -
                     loglik.plackett.tri(params=c(p[i,],params), counts[i,])
      }

      #control log(u) < log(r) (u in unif(0,1))
      if(log(runif(1)) < lar) {
        params[j] <- proposal[j]
        accept[it,3*M+j] <- 1
      } else {
        accept[it,3*M+j] <- 0
      }
      # --- Robbins-Monro update ---
      sigma[j] <- sigma[j] * exp((mcmc.options$delta / it) * (min(1,exp(lar)) - mcmc.options$target.accept))
      sigma.trace[it,j] <- sigma[j]
    }
    #SAVE THE CHAINS FOR EACH DATASET
    chains[it,] <- c(p,params)
    for (i in 1:M) {
      p_any[it,i]<-1-prob.fun.tri(c(p[i,],params))["p000"]
    }
  }
  if (mcmc.options$delta.lambda==0) {lambda.trace<-NULL}

  return(list(chains = chains,
              p_any = p_any,
              acceptance = colMeans(accept),
              sigma.trace = sigma.trace,
              lambda.trace = lambda.trace,
              counts = counts,
              priors = priors))
}
