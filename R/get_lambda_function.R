#' Construct spatio-temporal intensity functions with controlled separability
#'
#' Returns an intensity function \eqn{\lambda(x,y,t)} corresponding to one of four
#' models used for simulation experiments in Ghorbani et al. (2021). Each model is a
#' mixture of a separable "background" component and a structured (generally non-separable)
#' spatio-temporal Gaussian bump. The models provide different degrees of
#' space–time separability, allowing for controlled experiments on separability testing.
#'
#' The returned function is intended for use in simulation (e.g., for generating
#' spatio-temporal Poisson point patterns under varying degrees of separability).
#'
#' @param N Numeric scalar (> 0). Baseline intensity level (interpreted as an expected total count
#'   after scaling in the calling simulator; see details below).
#' @param g Numeric scalar (>= 0). Weight of the structured (non-separable) component.
#' @param model Integer in \code{1:4} indicating the  structure of the intensity function.
#' @param mu1 Numeric scalar. Mean of the temporal Gaussian background term (models 2 and 4).
#' @param sd1 Numeric scalar (> 0). Standard deviation of the temporal Gaussian background term (models 2 and 4).
#' @param mu2 Numeric vector of length 2. Mean of the spatial 2D Gaussian background term (models 3 and 4).
#' @param sd2 Numeric vector of length 2 with positive entries. Standard deviations of the spatial 2D Gaussian
#'   background term (models 3 and 4).
#' @param mu3 Numeric vector of length 3. Mean of the structured (non-separable) 3D Gaussian component.
#' @param sd3 Numeric vector of length 3 with positive entries. Standard deviations of the structured 3D Gaussian
#'   component.
#'
#' @details
#' The intensity is constructed as:
#' \deqn{\lambda(x,y,t) = \lambda_{\mathrm{bg}}(x,y,t) + g\, f_{st}(x,y,t),}
#' where \eqn{f_{st}} is a nonnegative 3D Gaussian density (via \code{\link{norm3d}}) and the background term
#' \eqn{\lambda_{\mathrm{bg}}} depends on \code{model}:
#' \describe{
#'   \item{1}{Homogeneous background: \eqn{\lambda_{\mathrm{bg}}(x,y,t) = (N-g)}}
#'   \item{2}{Temporal inhomogeneity only: \eqn{\lambda_{\mathrm{bg}}(x,y,t) = (N-g)\, f_t(t)}, where
#'            \eqn{f_t} is a 1D Gaussian density (\code{\link[stats]{dnorm}}).}
#'   \item{3}{Spatial inhomogeneity only: \eqn{\lambda_{\mathrm{bg}}(x,y,t) = (N-g)\, f_s(x,y)}, where
#'            \eqn{f_s} is a 2D Gaussian density (\code{\link{norm2d}}).}
#'   \item{4}{Separable spatial-temporal inhomogeneity:
#'            \eqn{\lambda_{\mathrm{bg}}(x,y,t) = (N-g)\, f_s(x,y)\, f_t(t)}.}
#' }
#'
#' Note: since Gaussian densities can exceed 1 for small standard deviations, \code{N} is best interpreted
#' as a scaling parameter used by the calling simulator. Ensure \eqn{\lambda(x,y,t)} is nonnegative over the
#' intended domain.
#'
#' See more details in Ghorbani et al. (2021), Section 6.1.
#'
#' @return A function of the form \code{function(x, y, t)} representing the selected intensity surface.
#'
#' @note
#' This function is primarily intended for generating intensity functions used in
#' simulation studies. In particular, \code{\link{rstpoispp}} calls
#' \code{get.lambda.function()} internally to construct intensity models for
#' simulating spatio-temporal Poisson point processes with controlled separability.
#'
#' @importFrom stats dnorm
#' @seealso
#' \code{\link{rstpoispp}},
#' \code{\link{norm3d}},
#' \code{\link{norm2d}},
#' \code{\link{chi2.test}},
#' \code{\link{dHS.test}}
#'
#' @author
#' Mohammad Ghorbani \email{mohammad.ghorbani@slu.se}\cr
#' Nafiseh Vafaei \email{nafiseh.vafaei@slu.se}
#'
#' @references
#' Ghorbani, M., Vafaei, N., Dvořák, J., and Myllymäki, M. (2021).
#' Testing the first-order separability hypothesis for spatio-temporal point patterns.
#' \emph{Computational Statistics & Data Analysis}, \bold{161}, 107245.
#'
#' @examples
#'
#' \donttest{
#' # Choose model 4: non-separable spatio-temporal intensity
#' lambda <- get.lambda.function(N = 210, g = 50, model = 4)
#' lambda(0.5, 0.5, 0.5)  # Evaluate intensity at center of space-time domain
#'
#' # Visualize spatial intensity at fixed time for model 2
#' lambda2 <- get.lambda.function(N = 200, g = 50, model = 2)
#' x <- y <- seq(0, 1, length.out = 100)
#' z <- outer(x, y, function(x, y) lambda2(x, y, t = 0.5))
#' par(mar = c(5, 4, 4, 6))
#' fields::image.plot(x, y, z, main = "Intensity at t = 0.5 (Model 2)", col = topo.colors(50))
#' }
#'
#' @export
get.lambda.function <- function(N, g, model = 1L,
                                mu1 = 0.5, sd1 = 0.2,
                                mu2 = c(0.5, 0.5), sd2 = c(0.2, 0.2),
                                mu3 = c(0.3, 0.3, 0.2), sd3 = c(0.05, 0.05, 0.05)) {

  # ---- validate scalars ----
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0) {
    stop("`N` must be a single finite positive numeric value.")
  }
  if (!is.numeric(g) || length(g) != 1L || !is.finite(g) || g < 0) {
    stop("`g` must be a single finite numeric value >= 0.")
  }
  model <- as.integer(model)
  if (!(model %in% 1:4)) stop("`model` must be an integer in 1:4.")

  # Keep your original constraint, but do not oversell what it guarantees.
  if (g > N) stop("`g` must be <= `N` to keep the background weight (N-g) nonnegative.")

  # ---- validate Gaussian parameters ----
  if (!is.numeric(mu1) || length(mu1) != 1L || !is.finite(mu1)) stop("`mu1` must be a finite numeric scalar.")
  if (!is.numeric(sd1) || length(sd1) != 1L || !is.finite(sd1) || sd1 <= 0) stop("`sd1` must be a finite positive numeric scalar.")

  if (!is.numeric(mu2) || length(mu2) != 2L || any(!is.finite(mu2))) stop("`mu2` must be a finite numeric vector of length 2.")
  if (!is.numeric(sd2) || length(sd2) != 2L || any(!is.finite(sd2)) || any(sd2 <= 0)) stop("`sd2` must be a finite positive numeric vector of length 2.")

  if (!is.numeric(mu3) || length(mu3) != 3L || any(!is.finite(mu3))) stop("`mu3` must be a finite numeric vector of length 3.")
  if (!is.numeric(sd3) || length(sd3) != 3L || any(!is.finite(sd3)) || any(sd3 <= 0)) stop("`sd3` must be a finite positive numeric vector of length 3.")

  # ---- define components ----
  bg_const <- N - g

  f_t <- function(t) stats::dnorm(t, mean = mu1, sd = sd1)
  f_s <- function(x, y) norm2d(x, y, mu = mu2, sd = sd2)
  f_st <- function(x, y, t) norm3d(x, y, t, mu = mu3, sd = sd3)

  # ---- build model functions ----
  funs <- list(
    # 1: homogeneous background
    function(x, y, t) bg_const + g * f_st(x, y, t),

    # 2: temporal background only
    function(x, y, t) bg_const * f_t(t) + g * f_st(x, y, t),

    # 3: spatial background only
    function(x, y, t) bg_const * f_s(x, y) + g * f_st(x, y, t),

    # 4: separable spatial-temporal background
    function(x, y, t) bg_const * f_s(x, y) * f_t(t) + g * f_st(x, y, t)
  )

  funs[[model]]
}
