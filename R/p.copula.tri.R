#' Trivariate Plackett Copula Cumulative Distribution Function
#'
#' Computes the value of a trivariate Plackett copula given
#' marginal values, bivariate dependence parameters, and a trivariate
#' dependence parameter.
#'
#' @param params A named numeric vector containing the arguments:
#'   \describe{
#'     \item{u}{Value in [0,1], first marginal prevalence}
#'     \item{v}{Value in [0,1], second marginal prevalence}
#'     \item{w}{Value in [0,1], third marginal prevalence}
#'     \item{phi.UV}{Dependence parameter for the (u,v) copula}
#'     \item{phi.UW}{Dependence parameter for the (u,w) copula}
#'     \item{phi.VW}{Dependence parameter for the (v,w) copula}
#'     \item{psi}{Trivariate dependence parameter}
#'   }
#'
#' @returns A numeric value in [0,1] corresponding to the trivariate copula
#'   cumulative distribution function evaluated at \eqn{(u,v,w)}.
#'   Returns `NA` if no valid solution exists or if parameters are invalid.
#'
#' @examples
#' # Example with arbitrary parameters
#' params <- c(u = 0.5, v = 0.6, w = 0.4,
#'               phi.UV = 2, phi.UW = 1.5, phi.VW = 0.8, psi=1.2)
#' p.copula.tri(params)
#' # Trivariate independence case
#' prob.fun.tri(c(params[-7],psi=1.2))
#'
#' @export
p.copula.tri <- function(params) {

  u <- params["u"]
  v <- params["v"]
  w <- params["w"]
  phi.UV <- params["phi.UV"]
  phi.UW <- params["phi.UW"]
  phi.VW <- params["phi.VW"]
  psi <- params["psi"]
  # Compute the bivariate copulas
  Cuv <- p.copula.bi(u, v, phi.UV)
  Cuw <- p.copula.bi(u, w, phi.UW)
  Cvw <- p.copula.bi(v, w, phi.VW)

  # a_i coefficients
  a1 <- Cvw
  a2 <- Cuw
  a3 <- Cuv
  a0 <- 1 - u - v - w + Cuv + Cuw + Cvw

  # b_i coefficients
  b1 <- Cuw + Cvw - w
  b2 <- Cuv + Cvw - v
  b3 <- Cuw + Cuv - u

  # Bounds for a valid solution
  vals.a <- c(a0, a1, a2, a3)
  vals.b <- c(0, b1, b2, b3)
  if (any(!is.finite(vals.a)) || any(!is.finite(vals.b))) {
    cat("Bound error",vals.a,vals.b,"\n")
    return(NA)
  }

  b <- max(vals.b)
  a <- min(vals.a)

  if (b > a) { ## no valid solution
    #cat("No valid solutions:",a,b,"\n")
    return(NA)
  } else if (b==a) { ## one solution?
    #cat("One solution:",a,"\n")
    return(a)
  }

  # Define the polynomial
  pol <- function(z) {
    psi*(a0 - z) * (a1 - z) * (a2 - z) * (a3 - z) -
      z * (z - b1) * (z - b2) * (z - b3)
  }
  # Find the root numerically
  sol <- uniroot(pol, interval = c(b, a))$root
  if (sol<b || sol>a) {cat("Uniroot error.\n")}
  #if (psi==1) {cat("psi=1:",u*v*w,sol,"\n")}
  return(sol)
}
