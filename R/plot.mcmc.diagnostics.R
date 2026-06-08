#' Plot MCMC Diagnostics for Plackett Copula Parameters
#'
#' This function makes diagnostic plots for MCMC runs of the Plackett copula parameters.
#' You can visualize the adaptation of proposal standard deviations (`sigma`),
#' trace plots of the chains (`trace`), or posterior histograms (`hist`) with MLE and prior overlays.
#' For details see \code{\link{mcmc.plackett}}.
#'
#' @param mcmc Output list from mcmc.plackett().
#' @param type Character. Type of plot to generate. Options are:
#'   \describe{
#'     \item{"trace"}{Trace plots of the MCMC chains for all parameters, with MLE overlaid.}
#'     \item{"hist"}{Posterior histograms of all parameters with MLE and prior distribution overlaid.}
#'     \item{"sigma"}{Trace plots of the proposal standard deviations for theta.xy, theta.xz, theta.yz, and psi.}
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
#' plot.mcmc.diagnostics(mcmc, type = "trace")
#'
plot.mcmc.diagnostics <- function(mcmc,type = c("trace", "hist", "sigma"),MLE=TRUE,log=TRUE) {
  type <- match.arg(type)

  priors<-mcmc$priors
  estimates <- get.estimates(mcmc$counts, MLE=MLE)

  #Easier access
  chains <- mcmc$chains
  if (log) {
    chains[,(ncol(chains)-3):ncol(chains)]<-log(chains[,(ncol(chains)-3):ncol(chains)])
    estimates[(ncol(chains)-3):ncol(chains)]<-log(estimates[(ncol(chains)-3):ncol(chains)])
    labs<-c(colnames(chains)[1:(ncol(chains)-4)],expression(log(phi[UV])),expression(log(phi[UW])),expression(log(phi[VW])),expression(log(psi)))
  } else {
    labs<-c(colnames(chains)[1:(ncol(chains)-4)],expression(phi[UV]),expression(phi[UW]),expression(phi[VW]),expression(psi))
  }

  leg<-rep("MLE",ncol(chains))
  if (!MLE) {
    leg[(ncol(chains)-3):ncol(chains)]<-"Robust estimate"
  }
  M<-nrow(mcmc$counts)

  if (type == "sigma") {
    labs<-c(expression(sigma[UV]),expression(sigma[UW]),expression(sigma[VW]),expression(sigma[psi]))
    for (j in 1:4) {
      plot(mcmc$sigma.trace[,j], type = "l", main = "Robbins-Monro scaling", ylab=labs[j],ylim=c(0,max(mcmc$sigma.trace[,j])),xlab="iteration")
    }

  } else if (type == "trace") {
    for (j in 1:ncol(chains)) {
      plot(chains[,j], type = "l", ylab = labs[j], main=paste("acceptance rate =",round(mcmc$acceptance[j],3)),xlab="iteration")
      abline(h = estimates[j], col = 2)
      if (j>ncol(chains)-4) {abline(h = 1-1*log, col="lightgrey")}
      legend("topright",leg[j],col = 2,lty=1)
    }

  } else if (type == "hist") {
    for (j in 1:ncol(chains)) {
      hist(chains[,j], breaks=40, prob=T, main = labs[j], col = "lightblue",xlab="")
      abline(v = estimates[j], col = 2)
      m0<-min(chains[,j])
      m1<-max(chains[,j])
      if (j<=M) {
        curve(dbeta(x, priors$a.u[j], priors$b.u[j]),from=m0,to=m1, add = TRUE, col = "darkgreen")
      } else if (j<=2*M) {
        curve(dbeta(x, priors$a.v[j-M], priors$b.v[j-M]),from=m0,to=m1, add = TRUE, col = "darkgreen")
      } else if (j<=3*M) {
        curve(dbeta(x, priors$a.w[j-2*M], priors$b.w[j-2*M]),from=m0,to=m1, add = TRUE, col = "darkgreen")
      } else if (!log) {
        curve(df(x, priors$d[j-3*M], priors$d[j-3*M]),from=m0,to=m1, add = TRUE, col = "darkgreen")
        abline(v = 1, col = "lightgrey")
      } else {
        abline(v = 0, col = "lightgrey")
      }
      legend("topright", c("Posterior", leg[j], "Prior"), fill=c("lightblue",0,0),border=c(1,0,0),lty = c(NA,1,1), col = c(NA,"red","darkgreen"),merge=TRUE)
    }
  }
}
