#' Plot Posteriors for Plackett Copula Parameters
#'
#' This function plots posteriors of the Plackett copula parameters.
#' For details see \code{\link{mcmc.plackett}}.
#'
#' @param mcmc Output list from mcmc.plackett().
#' @param type Character. Type of plot to generate. Options are:
#'   \describe{
#'     \item{"OR"}{Posteriors of the copula parameters as odds ratios.}
#'     \item{"prevalence"}{Posterior intervals for the prevalences, sorted by sample size.}
#'   }
#'
#' @return No return value. Generates plots for MCMC diagnostics.
#'
#' @export
#'
#' @examples
#' counts <- c(10, 5, 8, 6, 7, 4, 9, 11)
#' names(counts) <- c("n000","n001","n010","n011","n100","n101","n110","n111")
#' mcmc <- mcmc.plackett(counts, N = 1000)
#' plot.posteriors(mcmc)
#'
plot.posteriors <- function(mcmc,burnin=0,legs=rep("topleft",ncol(chains)),cols=c("lightblue","darkorange"),...) {
  #Easier access
  priors<-mcmc$priors
  chains <- mcmc$chains[(burnin+1):nrow(mcmc$chains),(ncol(mcmc$chains)-3):ncol(mcmc$chains)]
  labs<-c("Ascasris-Hookworm interaction","Ascaris-Trichuris interaction","Hookworm-Trichuris interaction","Ascaris-Hookworm-Trichuris interaction")
  xlabs<-c(expression(phi[UV]),expression(phi[UW]),expression(phi[VW]),expression(psi))
  M<-nrow(mcmc$counts)
  m0<-0
  m1<-max(chains)
  for (j in 1:4) {
    hist(chains[,j], breaks=seq(m0,m1,length.out=81), prob=T, main = labs[j],xlab=xlabs[j],xlim=c(m0,m1),col=cols[1],...)
    abline(v = 1, col = "lightgrey")
    curve(df(x, priors$d[j], priors$d[j]),from=m0,to=m1, add = TRUE, col = cols[2])
    legend(legs[j], c("Posterior", "Prior"),bg="white",fill=c(cols[1],0),border=c(1,0),lty = c(NA,1), col = c(NA,cols[2]),merge=TRUE)
  }
}

