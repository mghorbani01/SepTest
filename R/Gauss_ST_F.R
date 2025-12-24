#' Simulate a spatio-temporal Gaussian random field on a regular grid
#'
#' Simulates a space--time Gaussian random field on a regular \eqn{(x,y,t)} grid.
#' The field is returned as a 3D array and can be used as a latent field for
#' log-Gaussian Cox process (LGCP) simulation.
#'
#' The simulated field is a weighted sum of three independent Gaussian components:
#' \deqn{Z(x,y,t) = \sigma_1 Z_s(x,y) + \sigma_2 Z_t(t) + \sigma_3 Z_{st}(x,y,t),}
#' where \eqn{Z_s} is a purely spatial field, \eqn{Z_t} is a purely temporal field,
#' and \eqn{Z_{st}} is a spatio-temporal field with separable exponential covariance
#' in space and time.
#'
#' The function uses \code{\link[MASS]{mvrnorm}} for multivariate normal simulation
#' and \code{\link[fields]{rdist}} to compute pairwise distances for covariance
#' matrix construction.
#'
#' @param xlim,ylim,tlim Numeric vectors of length 2 giving the ranges for the spatial and temporal
#'   axes. Defaults are \code{c(0,1)} for each.
#' @param par1 Numeric vector of length 2 giving the temporal covariance parameters
#'   \code{c(variance, scale)} for an exponential covariance \eqn{var * exp(-d/scale)}.
#' @param par2 Numeric vector of length 2 giving the spatial covariance parameters
#'   \code{c(variance, scale)} for an exponential covariance \eqn{var * exp(-d/scale)}.
#' @param sigmas Numeric vector of length 3 specifying the weights
#' \eqn{(\sigma_1,\sigma_2,\sigma_3)} for combining the spatial, temporal, and
#' spatio-temporal components of the field.
#' @param grid Integer vector of length 3 giving the number of grid points in the \eqn{x}, \eqn{y}, and \eqn{t}
#'   directions.
#'
#' @return A list with components:
#' \describe{
#'   \item{Z}{Numeric array of dimension \code{c(nx, ny, nt)} containing simulated field values.}
#'   \item{xcoord}{Numeric vector of length \code{nx} with x-grid coordinates.}
#'   \item{ycoord}{Numeric vector of length \code{ny} with y-grid coordinates.}
#'   \item{tcoord}{Numeric vector of length \code{nt} with time-grid coordinates.}
#' }
#'
#' @details
#' Spatial and temporal covariances are exponential. The spatio-temporal component uses a separable
#' covariance \eqn{C_{st}((u,t),(u',t')) = C_s(u,u') C_t(t,t')}. Simulation is performed via Cholesky
#' factors without constructing the full \eqn{(nx*ny*nt) \times (nx*ny*nt)} covariance matrix.
#'
#' @author
#' Mohammad Ghorbani \email{mohammad.ghorbani@slu.se}\cr
#' Nafiseh Vafaei \email{nafiseh.vafaei@slu.se}
#'
#' @references
#' Ghorbani M., Vafaei N., Dvořák J., Myllymäki M. (2021).
#' Testing the first-order separability hypothesis for spatio-temporal point patterns.
#' \emph{Computational Statistics and Data Analysis}, \bold{161}, 107245.
#'
#' @seealso \code{\link[MASS]{mvrnorm}}, \code{\link[fields]{rdist}}
#'
#' @examples
#'
#' \donttest{
#' if (requireNamespace("MASS", quietly = TRUE) && requireNamespace("fields", quietly = TRUE)) {
#'   set.seed(1)
#'   field <- Gauss.st.F(
#'     xlim = c(0, 1), ylim = c(0, 1), tlim = c(0, 1),
#'     par1 = c(1, 0.05), par2 = c(1, 0.06),
#'     sigmas = c(0.5, 0.5, 1),
#'     grid = c(15, 15, 10)
#'   )
#' # Inspect dimensions and visualize one time slice
#' dim(field$Z)
#' image(field$xcoord, field$ycoord, field$Z[, , 1],
#'       main = "Gaussian Random Field (t = 1)",
#'       col = RColorBrewer::brewer.pal(11, "Spectral"))
#'   }
#' }
#'
#' @importFrom fields rdist
#' @importFrom MASS mvrnorm
#' @export
Gauss.st.F <- function(xlim = c(0, 1),
                       ylim = c(0, 1),
                       tlim = c(0, 1),
                       par1 = c(1, 0.05),
                       par2 = c(1, 0.06),
                       sigmas = c(0.5, 0.5, 1),
                       grid = c(15L, 15L, 10L)) {

  # ---- validate inputs ----
  if (!is.numeric(xlim) || length(xlim) != 2L || anyNA(xlim) || !(xlim[1] < xlim[2])) stop("`xlim` must be c(min,max) with min < max.")
  if (!is.numeric(ylim) || length(ylim) != 2L || anyNA(ylim) || !(ylim[1] < ylim[2])) stop("`ylim` must be c(min,max) with min < max.")
  if (!is.numeric(tlim) || length(tlim) != 2L || anyNA(tlim) || !(tlim[1] < tlim[2])) stop("`tlim` must be c(min,max) with min < max.")

  if (!is.numeric(par1) || length(par1) != 2L || anyNA(par1) || any(par1 <= 0)) stop("`par1` must be c(variance, scale) with both > 0.")
  if (!is.numeric(par2) || length(par2) != 2L || anyNA(par2) || any(par2 <= 0)) stop("`par2` must be c(variance, scale) with both > 0.")
  if (!is.numeric(sigmas) || length(sigmas) != 3L || anyNA(sigmas)) stop("`sigmas` must be a numeric vector of length 3.")
  if (!is.numeric(grid) || length(grid) != 3L || anyNA(grid) || any(grid <= 1)) stop("`grid` must be a vector of length 3 with values > 1.")
  if (any(abs(grid - round(grid)) > 0)) stop("`grid` must contain integer values.")
  grid <- as.integer(grid)

  nx <- grid[1L]; ny <- grid[2L]; nt <- grid[3L]
  ns <- nx * ny


  # Guard against huge covariance matrices in space
  # cov_xy is ns x ns. This becomes very large quickly.
  if (ns > 5000L) {
    stop("Spatial grid too large (nx*ny > 5000). Reduce `grid` to avoid huge covariance matrices.")
  }

  if (is.null(xlim)) xlim   = c(0, 1)
  if (is.null(ylim)) ylim   = c(0, 1)
  if (is.null(tlim)) tlim   = c(0, 1)

  x <- seq(xlim[1], xlim[2], length.out = grid[1])
  y <- seq(ylim[1], ylim[2], length.out = grid[2])
  t <- seq(tlim[1], tlim[2], length.out = grid[3])


  nx <- length(x)
  ny <- length(y)
  nt <- length(t)

  # Spatial coordinates
  xy <- as.matrix(expand.grid(x = x, y = y))

  # Temporal covariance (1D exponential)
  Dt <- fields::rdist(matrix(t), matrix(t))
  cov_t <- par1[1] * exp(-Dt / par1[2])
  z_t <- MASS::mvrnorm(1, mu = rep(0, nt), Sigma = cov_t)

  # Spatial covariance (2D exponential)
  Dxy <- fields::rdist(xy)
  cov_xy <- par2[1] * exp(-Dxy / par2[2])
  z_u <- MASS::mvrnorm(1, mu = rep(0, nx * ny), Sigma = cov_xy)
  z_u <- matrix(z_u, nrow = nx, ncol = ny)

  # Spatio-temporal fields (simulate nt spatial fields independently)
  z_st <- array(0, dim = c(nx, ny, nt))
  for (k in 1:nt) {
    temp_z <- MASS::mvrnorm(1, mu = rep(0, nx * ny), Sigma = cov_xy)
    z_st[,,k] <- matrix(temp_z, nrow = nx, ncol = ny)
  }

  # Combine components
  Z <- array(0, dim = c(nx, ny, nt))
  for (i in 1:nx) {
    for (j in 1:ny) {
      for (k in 1:nt) {
        Z[i,j,k] <- sigmas[1]*z_u[i,j] + sigmas[2]*z_t[k] + sigmas[3]*z_st[i,j,k]
      }
    }
  }

  return(list(Z=Z,xcoord=x,ycoord=y,tcoord=t))
}
