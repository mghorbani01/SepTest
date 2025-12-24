#' Bivariate normal density with independent components
#'
#' Evaluates the density of a bivariate normal distribution with mean vector
#' \eqn{\mu = (\mu_x,\mu_y)} and diagonal covariance matrix (independent components).
#' The density is the product of the two univariate normal densities:
#' \deqn{f(x,y) = \phi(x;\mu_x,\sigma_x)\,\phi(y;\mu_y,\sigma_y).}
#'
#' @param x Numeric vector of x-coordinate(s).
#' @param y Numeric vector of y-coordinate(s).
#' @param mu Numeric vector of length 2 giving \code{c(mu_x, mu_y)}.
#' @param sd Numeric vector of length 2 giving positive standard deviations \code{c(sd_x, sd_y)}.
#' @param log Logical; if \code{TRUE}, return the log-density.
#'
#' @return Numeric vector of densities (or log-densities) with length determined by standard
#' recycling rules for \code{x} and \code{y}.
#'
#' @author
#' Mohammad Ghorbani \email{mohammad.ghorbani@slu.se}\cr
#' Nafiseh Vafaei \email{nafiseh.vafaei@slu.se}
#'
#' @examples
#'
#' # Evaluate the density at the peak
#' norm2d(0.5, 0.5, mu = c(0.5, 0.5), sd = c(0.2, 0.2))
#'
#' # Evaluate at multiple x values
#' norm2d(c(0.3, 0.7), 0.5, mu = c(0.5, 0.5), sd = c(0.2, 0.2))
#'
#' # Visualize on a grid
#' x <- y <- seq(0, 1, length.out = 100)
#' f <- Vectorize(function(x, y) norm2d(x, y, mu = c(0.5, 0.5), sd = c(0.2, 0.2)))
#' z <- outer(x, y, f)
#' image(x, y, z, col = terrain.colors(50), main = "Bivariate Normal Intensity")
#' contour(x, y, z, add = TRUE)
#'
#' @export
norm2d<- function(x, y, mu = c(0, 0), sd = c(1, 1), log = FALSE) {

  if (!is.numeric(x) || !is.numeric(y)) stop("`x` and `y` must be numeric.")
  if (!is.numeric(mu) || length(mu) != 2L || anyNA(mu) || any(!is.finite(mu))) {
    stop("`mu` must be a finite numeric vector of length 2.")
  }
  if (!is.numeric(sd) || length(sd) != 2L || anyNA(sd) || any(!is.finite(sd)) || any(sd <= 0)) {
    stop("`sd` must be a finite positive numeric vector of length 2.")
  }
  if (!is.logical(log) || length(log) != 1L) stop("`log` must be a single logical value.")

  lx <- stats::dnorm(x, mean = mu[1L], sd = sd[1L], log = log)
  ly <- stats::dnorm(y, mean = mu[2L], sd = sd[2L], log = log)

  if (log) lx + ly else lx * ly
}





