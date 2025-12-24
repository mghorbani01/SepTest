#' Upper bound for spatio-temporal intensity models
#'
#' Computes a practical upper bound for the spatio-temporal intensity models used in
#' \code{\link{get.lambda.function}}. The bound is intended for thinning/rejection sampling
#' in simulation routines such as \code{\link{rstpoispp}}.
#'
#' The bound is computed using analytic maxima of Gaussian density components (at their modes),
#' which yields a conservative and fast-to-evaluate upper bound when the component functions are
#' Gaussian product densities.
#'
#' @param N Numeric scalar (> 0). Total expected number of events (baseline intensity level).
#' @param g Numeric scalar (>= 0). Weight of the structured (non-separable) component.
#' @param model Integer in \code{1:4}. See \code{\link{get.lambda.function}}.
#' @param mu1 Numeric scalar. Mean of the temporal Gaussian background term (models 2 and 4).
#' @param sd1 Numeric scalar (> 0). Standard deviation of the temporal Gaussian background term.
#' @param mu2 Numeric vector of length 2. Mean of the spatial Gaussian background term (models 3 and 4).
#' @param sd2 Numeric vector of length 2 with positive entries. Standard deviations of the spatial background term.
#' @param mu3 Numeric vector of length 3. Mean of the structured spatio-temporal Gaussian component.
#' @param sd3 Numeric vector of length 3 with positive entries. Standard deviations of the structured component.
#'
#' @details
#' The intensity models are mixtures of a background term and a structured spatio-temporal Gaussian bump.
#' This function returns an upper bound obtained by evaluating each Gaussian density component at its mode.
#' This upper bound is typically sufficient for rejection sampling when generating
#' realizations of inhomogeneous Poisson or Cox point processes.
#'
#' If \code{norm2d} and \code{norm3d} in this package are Gaussian product densities (independent components),
#' then the maxima are available in closed form:
#' \itemize{
#'   \item \eqn{\max_t \phi(t;\mu,\sigma) = 1/(\sigma\sqrt{2\pi})}
#'   \item \eqn{\max_{x,y} \phi(x;\mu_x,\sigma_x)\phi(y;\mu_y,\sigma_y) = 1/(2\pi\sigma_x\sigma_y)}
#'   \item \eqn{\max_{x,y,t} \prod_{d=1}^3 \phi(\cdot;\mu_d,\sigma_d) = 1/((2\pi)^{3/2}\sigma_x\sigma_y\sigma_t)}
#' }
#'
#' If your \code{norm2d}/\code{norm3d} use a different parameterization, this bound should be updated accordingly.
#'
#' @return Numeric scalar. A conservative upper bound for the selected intensity model.
#'
#' @importFrom stats dnorm
#'
#' @references
#' Ghorbani, M., Vafaei, N., Dvořák, J., and Myllymäki, M. (2021).
#' Testing the first-order separability hypothesis for spatio-temporal point patterns.
#' \emph{Computational Statistics and Data Analysis}, \bold{161}, 107245.
#'
#' @seealso \code{\link{get.lambda.function}}, \code{\link{rstpoispp}}
#'
#' @author
#' Nafiseh Vafaei \email{nafiseh.vafaei@slu.se}\cr
#' Mohammad Ghorbani \email{mohammad.ghorbani@slu.se}
#'
#' @examples
#' # Example 1: Homogeneous model (Model 1)
#' get.lambda.max(N = 200, g = 50, model = 1)
#'
#' # Example 2: Non-separable spatio-temporal model (Model 4)
#' get.lambda.max(N = 200, g = 50, model = 4)
#'
#' @export
get.lambda.max <- function(N, g, model = 1L,
                           mu1 = 0.5, sd1 = 0.2,
                           mu2 = c(0.5, 0.5), sd2 = c(0.2, 0.2),
                           mu3 = c(0.3, 0.3, 0.2), sd3 = c(0.05, 0.05, 0.05)) {

  # ---- validate ----
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0) stop("`N` must be a finite positive numeric scalar.")
  if (!is.numeric(g) || length(g) != 1L || !is.finite(g) || g < 0) stop("`g` must be a finite numeric scalar >= 0.")
  model <- as.integer(model)
  if (!(model %in% 1:4)) stop("`model` must be an integer in 1:4.")
  if (g > N) stop("`g` must be <= `N` so that (N - g) is nonnegative.")

  if (!is.numeric(sd1) || length(sd1) != 1L || !is.finite(sd1) || sd1 <= 0) stop("`sd1` must be a finite positive numeric scalar.")
  if (!is.numeric(sd2) || length(sd2) != 2L || any(!is.finite(sd2)) || any(sd2 <= 0)) stop("`sd2` must be a finite positive numeric vector of length 2.")
  if (!is.numeric(sd3) || length(sd3) != 3L || any(!is.finite(sd3)) || any(sd3 <= 0)) stop("`sd3` must be a finite positive numeric vector of length 3.")

  # ---- analytic maxima of Gaussian product densities ----
  max_1d <- function(sd) 1 / (sd * sqrt(2 * pi))
  max_2d <- function(sd2) 1 / (2 * pi * sd2[1L] * sd2[2L])
  max_3d <- function(sd3) 1 / ((2 * pi)^(3/2) * sd3[1L] * sd3[2L] * sd3[3L])

  bg_const <- N - g
  max_st <- max_3d(sd3)

  # background maxima depending on model
  max_bg <- switch(
    as.character(model),
    "1" = bg_const,
    "2" = bg_const * max_1d(sd1),
    "3" = bg_const * max_2d(sd2),
    "4" = bg_const * max_2d(sd2) * max_1d(sd1)
  )

  # mixture upper bound
  max_bg + g * max_st
}
