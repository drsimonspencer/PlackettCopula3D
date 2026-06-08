#' CDF of the Bivariate Plackett Copula
#'
#' Computes the cumulative distribution function (CDF) of the bivariate
#' Plackett copula for given arguments \code{u}, \code{v} and dependence
#' parameter \code{phi}.
#'
#' @param u Numeric scalar in [0,1]. First marginal prevalence.
#' @param v Numeric scalar in [0,1]. Second marginal prevalence.
#' @param phi Positive numeric scalar. Dependence parameter of the
#' Plackett copula. Values close to 1 correspond to near-independence.
#'
#' @return A numeric scalar: the value of the copula CDF at (\code{u}, \code{v})
#' given \code{phi}. Returns \code{NaN} if inputs are invalid.
#'
#' @examples
#' # Independent case (phi = 1)
#' p.copula.bi(0.5, 0.5, phi = 1)
#'
#' # Moderate dependence
#' p.copula.bi(0.3, 0.7, phi = 2)
#'
#' # Edge cases
#' p.copula.bi(0, 0, phi = 1.5)
#' p.copula.bi(1, 1, phi = 0.5)
#'
#' @export
p.copula.bi <- function(u, v, phi) {


  if (!is.finite(phi) || phi <= 0 || !is.finite(u) || u < 0 || u > 1 ||
      !is.finite(v) || v < 0 || v > 1) {
    return(NaN)
  }

  if (abs(phi - 1) < 1e-8) {
    return(u * v)
  } else {
    num <- 1 + (phi - 1) * (u + v)
    delta <- num^2 - 4 * phi * (phi - 1) * u * v
    if (!is.finite(delta) || delta < 0) return(NaN)
    sqrt.term <- sqrt(delta)
    out <- (num - sqrt.term) / (2 * (phi - 1))
    names(out)<-NULL
    if (is.na(out) || out<max(0,u+v-1) || out>min(u,v)) {cat("Bivariate copula error:",u,v,phi,out,"\n")}
    return(out)
  }
}
