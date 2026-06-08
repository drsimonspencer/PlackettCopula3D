#' Joint probabilities for the trivariate Plackett copula
#'
#' Computes the eight joint probabilities \eqn{P(U=i, V=j, W=k)},
#' implied by the trivariate Plackett copula.
#' The function uses as inputs the marginal prevalences, the three
#' bivariate Plackett parameters, and the trivariate dependence parameter \code{psi}.
#'
#' @param params A named numeric vector with the parameters:
#'   \itemize{
#'     \item \code{u, v, w} : marginal prevalences for the three binary variables
#'     \item \code{phi.UV, phi.UW, phi.VW} : bivariate Plackett copula parameters
#'       for the corresponding pairs
#'     \item \code{psi} : dependence parameter for trivariate copula
#'   }
#'
#' @returns A named numeric vector of length 8 with entries:
#' \code{p000, p001, p010, p011, p100, p101, p110, p111}, corresponding to
#' the joint probabilities of the trivariate distribution.
#'
#' @examples
#' params <- c(u = 0.5, v = 0.5, w = 0.5,
#'               phi.UV = 1, phi.UW = 1, phi.VW = 1, psi=1)
#' prob.fun.tri(params)
#'
#' @export
prob.fun.tri <- function(params) {

  #compute the copulas (bivariate and trivariate)
  Cuv <- p.copula.bi(params["u"], params["v"], params["phi.UV"])
  Cuw <- p.copula.bi(params["u"], params["w"], params["phi.UW"])
  Cvw <- p.copula.bi(params["v"], params["w"], params["phi.VW"])
  Cuvw <- p.copula.tri(params)

  out<-c(1-params["u"]-params["v"]-params["w"] +Cuv +Cuw +Cvw-Cuvw,
         params["w"] - Cuw - Cvw + Cuvw,
         params["v"] - Cuv - Cvw + Cuvw,
         Cvw - Cuvw,
         params["u"] - Cuv - Cuw + Cuvw,
         Cuw - Cuvw,
         Cuv - Cuvw,
         Cuvw)
  names(out)<-c("p000",
                "p001",
                "p010",
                "p011",
                "p100",
                "p101",
                "p110",
                "p111")

  return(out)
}
