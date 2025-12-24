#' Simulate a Spatio-Temporal Determinantal Point Process (DPP) Based on Spectral Density
#'
#' Generates a realization of a spatio-temporal determinantal point process (DPP)
#' using a user-defined spectral density model.
#' The function supports both separable (\code{model = "S"}) and non-separable
#' (\code{model = "NS"}) dependence structures and allows for either exponential
#' or Mat\'ern-type spectral densities.
#' The simulation is performed on a 3D spatio-temporal frequency grid and can
#' incorporate user-specified intensity functions for thinning.
#'
#'@param mode Character. Type of dependence model:
#' \describe{
#'   \item{\code{"Stationary"}}{ A homogeneous spatio-temporal DPP is generated.}
#'   \item{\code{"Inhomogeneous"}}{A non-homogeneous spatio-temporal DPP is generated.}
#' }
#' @param model Character. Type of dependence model:
#' \describe{
#'   \item{\code{"S"}}{Separable spatio-temporal covarince function model, where space and time components are separable.}
#'   \item{\code{"NS"}}{Non-separable spatio-temporal model, allowing interaction between space and time.}
#' }
#' @param spectral Character. Type of spectral density function to use:
#' \code{"exp"} for exponential spectral form or \code{"matern"} for a Mat\'ern-type spectral model.
#' @param alpha_s Numeric. Spatial decay or range parameter in the spectral density.
#' @param alpha_t Numeric. Temporal decay or range parameter in the spectral density.
#' @param nu Numeric. Smoothness parameter for the Mat\'ern-type spectral density
#' (only relevant if \code{spectral = "matern"}). Default is \code{2}.
#' @param eps Numeric. Degree of separability for the Mat\'ern model:
#' \code{eps = 0} corresponds to full non-separability,
#' \code{eps = 1} yields complete separability,
#' and intermediate values provide partial separability.
#' @param lambda_s Optional. Intensity function \code{lambda_s(u, t)} for separable models.
#' If not provided, a default function is used.
#' @param lambda_non_s Optional. Intensity function \code{lambda_non_s(u, t)} for
#' non-separable models. If not provided, a default function is used.
#' @param lambda_max Numeric. The maximum intensity value used for thinning.
#' Must be specified.
#' @param grid_size Numeric. Half-width of the spatio-temporal frequency grid.
#' The total grid size is \code{(2 * grid_size + 1)^3}.
#'
#' @details
#' This function implements a spectral simulation method for spatio-temporal DPPs,
#' following the theoretical framework introduced in
#' Vafaei et al. (2023) and  Ghorbani et al. (2025).
#'
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Construct a 3D grid of spatial and temporal frequency components
#'   \eqn{(\omega_x, \omega_y, \tau)}.
#'   \item Evaluate the chosen spectral density \eqn{\phi(\omega, \tau)} across the grid.
#'   \item Use the resulting spectral values as eigenvalues to simulate a realization
#'   of a DPP via \eqn{spatstat.model::rdpp()}.
#'   \item Optionally apply thinning using a user-defined intensity function
#'   \eqn{\lambda(u, t)}, scaled by \code{lambda_max}, to induce inhomogeneity.
#' }
#'
#' Two spectral families are supported:
#' \itemize{
#'   \item \strong{Exponential form:}
#'     \deqn{
#'       \phi(\omega, \tau) \propto
#'         \exp\left[-(\pi \alpha_s |\omega|)^2\right]
#'         \left(1 + 4 (\pi \alpha_t \tau)^2\right)^{-1}.
#'     }
#'   \item \strong{Mat\'ern-type form:}
#'     \deqn{
#'       \phi_{\epsilon}(\omega, \tau) \propto
#'         \left(
#'           \alpha_s^2 \alpha_t^2
#'           + \alpha_t^2 |\omega|^2
#'           + \alpha_s^2 \tau^2
#'           + \epsilon |\omega|^2 \tau^2
#'         \right)^{-\nu},
#'     }
#'     where \eqn{\epsilon \in [0, 1]} determines the degree of separability between
#'     space and time.
#' }
#'
#' This framework enables simulation of spatio-temporal point patterns that exhibit
#' varying degrees of spatial–temporal dependence, providing a versatile tool for
#' evaluating separability tests and modeling non-separable dynamics.
#'
#' @return
#' A numeric matrix with three columns (\code{x}, \code{y}, \code{t}) representing
#' the retained spatio-temporal events after thinning.
#'
#' @seealso
#' \code{\link{plot_stpp}} for visualizing spatio-temporal point patterns.
#'
#' @references
#' Vafaei, N., Ghorbani, M., Ganji, M., and Myllymäki, M. (2023).
#' Spatio-temporal determinantal point processes.
#' *arXiv:2301.02353*.
#'
#' Ghorbani, M., Vafaei, N., and Myllymäki, M. (2025).
#' A kernel-based test for the first-order separability of spatio-temporal point processes.
#' \emph{TEST}, 34, 580-611.
#' https://doi.org/10.1007/s11749-025-00972-y
#'
#' @examples
#'
#' \donttest{
#'
#' # Simulate a stationary separable Mat\'ern ST-DPP
#'  if (requireNamespace("spatstat", quietly = TRUE)) {
#' sim <- rstDPP(
#'   mode     = "stationary",
#'   model    = "S",
#'   spectral = "matern",
#'   alpha_s  = 10,
#'   alpha_t  = 4.7,
#'   nu       = 2,
#'   eps      = 1,
#'   lambda_max = 70,
#'   grid_size  = 2
#' )
#'  plot_stDPP(sim, type = "3D",  alpha_s = 10, alpha_t = 4.7)
#' # example 2
#' # Generate realization
#' sim <- rstDPP(mode = "stationary",
#'                 model = "S",
#'                 spectral = "matern",
#'                 alpha_s = 10, alpha_t = 4.7,
#'                 nu = 2,
#'                eps = 1,
#'                lambda_s = 70,
#'                lambda_non_s = NULL,
#'                grid_size = 2,
#'                lambda_max=70)
#' head(sim)
#' }
#'}
#'
#' @importFrom stats runif na.omit
#' @importFrom spatstat.model rdpp
#' @export
rstDPP <- function(mode = c("stationary","inhom"),
                            model = c("S","NS"),
                            spectral = c("exp","matern"),
                            alpha_s, alpha_t=NULL,
                            lambda_max,
                            nu = 2, eps=1,
                            lambda_s = NULL,
                            lambda_non_s = NULL,
                            grid_size = 4) {

  if (!requireNamespace("spatstat", quietly = TRUE)) stop("Package 'spatstat' is required.")
  mode     <- match.arg(mode)
  model    <- match.arg(model)
  spectral <- match.arg(spectral)

  if (missing(alpha_s) || missing(lambda_max))
    stop("Please provide BOTH alpha_s and lambda_max (numeric).")

  # Define rho_max once; used by both default intensities
  rho_max <- lambda_max

  # Default inhomogeneous intensities for thinning
  if (mode == "inhom" && model == "S" && is.null(lambda_s)) {
    lambda_s <- function(x, y, t) rho_max * pmax(0, t) * (1 + cos(x + y))
  }
  if (mode == "inhom" && model == "NS" && is.null(lambda_non_s)) {
    lambda_non_s <- function(x, y, t) rho_max * (1 + cos(5 * (x + y + t)))
  }

  # Default stationary intensities
  if (mode == "stationary" && model == "S" ) {
    lambda_s <- function(x, y, t) rho_max
  }
  if (mode == "stationary" && model == "NS" ) {
    lambda_non_s <- function(x, y, t) rho_max
  }

  # -------- derive alpha_t from (alpha_s, lambda_max) only if missing --------
  if (spectral == "exp") {
    if (missing(alpha_t) || is.null(alpha_t)) {
      alpha_t <- 1 / (2 * pi * lambda_max * alpha_s^2)
    }
  } else { # "matern"
    if (missing(alpha_t) || is.null(alpha_t)) {
      if (eps == 1) {
        alpha_t <- (4 * lambda_max) / (pi^2 * alpha_s^2)
      } else if (eps == 0) {
        alpha_t <- (2 * lambda_max) / (pi^2 * alpha_s^2)
      } else {
        stop("eps must be 0 (non-separable) or 1 (separable) for the Matern case.")
      }
    }
  }

  if (!is.numeric(alpha_t) || alpha_t <= 0) stop("alpha_t must be a positive number.")

  # -------- spectral grid --------
  a <- grid_size
  Index <- expand.grid(x = -a:a, y = -a:a, t = -a:a)
  norm.Omeg <- sqrt(Index$x^2 + Index$y^2)

  # --- spectral densities (do not overwrite lambda_max) ---
  # if (spectral == "exp") {
  #   specden.st <- function(x, t, lambda, alpha_s, alpha_t) {
  #     phi_s <- pi * alpha_s^2 * exp(-(pi * alpha_s * x)^2)
  #     phi_t <- (2 * alpha_t) / (1 + 4 * (pi * alpha_t * t)^2)
  #     lambda * phi_s * phi_t
  #   }
  #   eig <- specden.st(norm.Omeg, Index$t, lambda = lambda_max, alpha_s = alpha_s, alpha_t = alpha_t)
  #
  # } else { # Matern (nu used; eps controls separability)
  #   specden.st <- function(alpha_s, alpha_t, eps, nu, x, t) {
  #     g <- alpha_s^4 * alpha_t^4
  #     g * (alpha_s^2 * alpha_t^2 + alpha_t^2 * x^2 + alpha_s^2 * t^2 + eps * x^2 * t^2)^(-nu)
  #   }
  #   eig <- specden.st(alpha_s, alpha_t, eps, nu, norm.Omeg, Index$t)
  # }


  # -------- unified spectral density function --------
  specden.st <- function(x, t, alpha_s, alpha_t, spectral = c("exp","matern"), nu = 2, eps = 1, lambda = 1) {
    spectral <- match.arg(spectral)
    if (spectral == "exp") {
      phi_s <- pi * alpha_s^2 * exp(-(pi * alpha_s * x)^2)
      phi_t <- (2 * alpha_t) / (1 + 4 * (pi * alpha_t * t)^2)
      return(lambda * phi_s * phi_t)
    } else { # matern
      g <- alpha_s^4 * alpha_t^4
      return(g * (alpha_s^2 * alpha_t^2 + alpha_t^2 * x^2 + alpha_s^2 * t^2 + eps * x^2 * t^2)^(-nu))
    }
  }

  eig <- specden.st(
    x = norm.Omeg,
    t = Index$t,
    alpha_s = alpha_s,
    alpha_t = alpha_t,
    spectral = spectral,
    nu = nu,
    eps = eps,
    lambda = lambda_max
  )
  if (anyNA(eig)) warning("NA values in eigenvalues were removed.")
  eig <- as.vector(na.omit(eig))

  # -------- simulate from DPP --------
  X <- spatstat.model::rdpp(eig, Index)  # assumes eigenvalues are within [0,1] as required by rdpp()
  x_val <- X$data$x
  y_val <- X$data$y
  t_val <- X$data$z  # ensure your rdpp() returns 'z' for time; otherwise use X$data$t

  # -------- inhomogeneous thinning --------
  intensity <- if (model == "S") lambda_s(x_val, y_val, t_val) else lambda_non_s(x_val, y_val, t_val)
  p <- pmin(1, pmax(0, intensity / lambda_max))
  retain <- runif(length(p)) <= p

  xyt <- cbind(x_val, y_val, t_val)
  X_final <- xyt[retain, , drop = FALSE]

  # attributes
  attr(X_final, "alpha_s")    <- alpha_s
  attr(X_final, "alpha_t")    <- alpha_t
  attr(X_final, "lambda_max") <- lambda_max
  attr(X_final, "spectral")   <- spectral
  attr(X_final, "model")      <- model
  attr(X_final, "eps")        <- eps
  attr(X_final, "nu")         <- nu

  return(X_final)
}





