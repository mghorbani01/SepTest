#' Simulate a Spatio-Temporal Log-Gaussian Cox Process (LGCP)
#'
#' Generates a realization of a spatio-temporal LGCP over a user-defined domain. The process is simulated
#' using a log-Gaussian random field combined with a deterministic trend function, and points are generated
#' by thinning a homogeneous Poisson process.
#'
#' @param xlim,ylim,tlim Numeric vectors of length 2 specifying the spatial and temporal domains.
#' @param grid Integer vector of length 3 specifying the number of grid cells in x, y, and t.
#' @param mu Optional. A function of (x, y, t, par) defining a deterministic trend. Default is nonlinear.
#' @param Lambda Optional. A user-supplied 3D intensity array or function. If \code{NULL}, it's generated from the latent Gaussian field.
#' @param Lmax Optional. Maximum intensity used for thinning. Can be numeric or a function. If \code{NULL}, it's computed automatically.
#' @param par1,par2 Parameters for temporal and spatial exponential covariance models, respectively.
#' @param sigmas Weights for combining spatial, temporal, and spatio-temporal components of the latent Gaussian field.
#' @param mu_par Parameters passed to the default trend function \code{mu()} if not user-supplied.
#'
#' @return A list with:
#' \describe{
#'   \item{st.lgcp}{A data frame of simulated spatio-temporal points.}
#'   \item{RF}{The latent Gaussian field output from \code{\link{Gauss.st.F}}.}
#' }
#' @importFrom stats rpois runif
#'
#' @examples
#'\donttest{
#' out <- rstLGCPP(xlim = c(0,1),
#'                 ylim = c(0,1),
#'                 tlim = c(0,1),
#'                 grid = c(15,15,10))
#' plot_stlgcp(data = out)
#' plot_stpp(data = out$st.lgcp, type = "3D")
#' }
#' @export

rstLGCPP <- function(xlim=NULL,ylim=NULL,tlim=NULL,grid= c(15,15,10),
                     mu=NULL,
                     Lambda = NULL,
                     Lmax = NULL,
                     par1 = c(1, 0.05),
                     par2 = c(1, 0.06),
                     sigmas = c(0.5, 0.5, 1),
                     mu_par = c(1.2, 0.25, 5)) {

  # Deterministic trend function
  if (is.null(mu))  {mu <- function(x, y, t, par) {
    par[1] + par[2] * (x - t) + par[3] * x * t
  }
  }
  # simulate or accept a user‐supplied intensity array Lambda_array

  # Gaussian random field on the grid
  g <- Gauss.st.F(xlim = xlim, ylim = ylim, tlim = tlim,
                  par1 = par1, par2 = par2, sigmas = sigmas,grid=grid)

  x=g$xcoord; y=g$ycoord; t=g$tcoord
  # deterministic trend mu(x,y,t)
  nx <- length(x); ny <- length(y); nt <- length(t)
  mu_arr <- array(NA, dim = c(nx, ny, nt))
  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      for (k in seq_along(t)) {
        mu_arr[i, j, k] <- mu(x[i], y[j], t[k], par = mu_par)
      }
    }
  }

  if (is.null(Lambda)) {
    Lambda_array <- exp(g$Z + mu_arr)
  } else if (is.function(Lambda)) {
    # user‐supplied function to compute intensity
    Lambda_array <- Lambda(x, y, t)
  } else {
    # user‐supplied intensity array
    Lambda_array <- Lambda
  }

  # determine the upper bound for thinning
  if (is.null(Lmax)) {
    Lmax_val <- max(Lambda_array) +
      0.05 * (max(Lambda_array) - min(Lambda_array))
  } else if (is.function(Lmax)) {
    # user‐supplied function to compute Lmax from the array
    Lmax_val <- Lmax(Lambda_array)
  } else {
    Lmax_val <- Lmax
  }

  n <- rpois(1, Lmax_val)
  XP <- data.frame(x = runif(n), y = runif(n), t = runif(n))
  X.P <- XP[order(XP$t), ]

  xstep <- x[2] - x[1]
  ystep <- y[2] - y[1]
  tstep <- t[2] - t[1]

  lambdaX <- numeric(0)

  for (k in 1:n) {
    indX <- (x - xstep/2 <= X.P$x[k]) & (X.P$x[k] < x + xstep/2)
    indY <- (y - ystep/2 <= X.P$y[k]) & (X.P$y[k] < y + ystep/2)
    indT <- (t - tstep/2 <= X.P$t[k]) & (X.P$t[k] < t + tstep/2)

    ix <- which(indX)
    iy <- which(indY)
    it <- which(indT)

    if (length(ix) == 1 && length(iy) == 1 && length(it) == 1) {
      lambdaX <- c(lambdaX, Lambda_array[ix, iy, it])
    } else {
      lambdaX <- c(lambdaX, 0)
    }
  }

  prob <- lambdaX / Lmax_val
  retain <- runif(n) <= prob
  st.lgcp <- data.frame(x = X.P$x[retain], y = X.P$y[retain], t = X.P$t[retain])
  return(list(st.lgcp=st.lgcp, RF=g))
}
