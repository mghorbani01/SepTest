#' Trivariate normal density with independent components in space-time
#'
#' Evaluates a trivariate normal density on \eqn{(x,y,t)} with independent components
#' (diagonal covariance). The density is the product of three univariate normal densities:
#' \deqn{f(x,y,t) = \phi(x;\mu_x,\sigma_x)\,\phi(y;\mu_y,\sigma_y)\,\phi(t;\mu_t,\sigma_t).}
#'
#' @param x Numeric vector of x-coordinate(s).
#' @param y Numeric vector of y-coordinate(s).
#' @param t Numeric vector of time coordinate(s).
#' @param mu Numeric vector of length 3 giving \code{c(mu_x, mu_y, mu_t)}.
#' @param sd Numeric vector of length 3 giving positive standard deviations \code{c(sd_x, sd_y, sd_t)}.
#' @param log Logical; if \code{TRUE}, return the log-density.
#'
#' @return Numeric vector of densities (or log-densities) with length determined by standard
#' recycling rules for \code{x}, \code{y}, and \code{t}.
#'
#' @seealso \code{\link{norm2d}}, \code{\link{get.lambda.function}}, \code{\link{st.intensity}}
#'
#' @author
#' Mohammad Ghorbani \email{mohammad.ghorbani@slu.se}\cr
#'
#' @references
#' Ghorbani, M., Vafaei, N., Dvořák, J., and Myllymäki, M. (2021).
#' Testing the first-order separability hypothesis for spatio-temporal point patterns.
#' \emph{Computational Statistics & Data Analysis}, \bold{161}, 107245.
#'
#' @examples
#' norm3d(0.3, 0.3, 0.2)  # peak value at the mean (with default parameters)
#' norm3d(c(0.2, 0.3), 0.3, 0.2)
#'
#' \donttest{
#' x <- y <- seq(0, 1, length.out = 100)
#' z <- outer(x, y, function(x, y) norm3d(x, y, t = 0.2))
#' image(x, y, z, col = heat.colors(50), main = "Spatial slice of norm3d at t = 0.2")
#' }
#'
#' @export
norm3d <- function(x, y, t,
                   mu = c(0.3, 0.3, 0.2),
                   sd = c(0.05, 0.05, 0.05),
                   log = FALSE) {

  if (!is.numeric(x) || !is.numeric(y) || !is.numeric(t)) {
    stop("`x`, `y`, and `t` must be numeric.")
  }
  if (!is.numeric(mu) || length(mu) != 3L || anyNA(mu) || any(!is.finite(mu))) {
    stop("`mu` must be a finite numeric vector of length 3.")
  }
  if (!is.numeric(sd) || length(sd) != 3L || anyNA(sd) || any(!is.finite(sd)) || any(sd <= 0)) {
    stop("`sd` must be a finite positive numeric vector of length 3.")
  }
  if (!is.logical(log) || length(log) != 1L) stop("`log` must be a single logical value.")

  lx <- stats::dnorm(x, mean = mu[1L], sd = sd[1L], log = log)
  ly <- stats::dnorm(y, mean = mu[2L], sd = sd[2L], log = log)
  lt <- stats::dnorm(t, mean = mu[3L], sd = sd[3L], log = log)

  if (log) lx + ly + lt else lx * ly * lt
}
