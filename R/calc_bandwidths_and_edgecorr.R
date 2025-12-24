#' Compute bandwidths and (Gaussian) edge-correction factors for spatio-temporal kernel intensity estimation
#'
#' Computes spatial and temporal bandwidths for kernel-based estimation of the intensity of a
#' spatio-temporal point pattern. The spatial bandwidth is estimated from the spatial coordinates
#' using Diggle's (1985) mean-square error method via \code{\link[spatstat.explore]{bw.diggle}}.
#' The temporal bandwidth is estimated from the time coordinates using the Sheather--Jones
#' direct plug-in method (\code{bw = "SJ-dpi"}) as implemented in \code{\link[stats]{density}}.
#'
#' Edge-correction factors are computed for Gaussian kernels as the kernel mass inside the
#' observation window: \eqn{c(x) = \int_W K_h(u-x)\,du}. The temporal correction is computed
#' exactly on \code{t.region}. The spatial correction is computed using a **bounding-box
#' approximation** of the polygonal spatial window (i.e., \code{W} is replaced by its bounding
#' rectangle).
#'
#' @param X A numeric matrix or data frame with at least three columns giving the \eqn{x}-, \eqn{y}-,
#'   and time coordinates \eqn{t} of observed spatio-temporal events. Each row corresponds to one event.
#' @param s.region A numeric matrix with two columns giving the vertices of the polygonal spatial
#'   observation window in order (the first vertex need not be repeated at the end).
#' @param t.region A numeric vector of length 2 giving the temporal observation window
#'   \code{c(tmin, tmax)} with \code{tmin < tmax}.
#' @param n.grid An integer vector of length 3 giving the number of grid cells in the \eqn{x}, \eqn{y},
#'   and time dimensions. Only \code{n.grid[3]} is used when estimating the temporal bandwidth.
#' @param epsilon Optional positive numeric. Spatial bandwidth. If \code{NULL}, estimated using
#'   \code{bw.diggle}.
#' @param delta Optional positive numeric. Temporal bandwidth. If \code{NULL}, estimated using
#'   \code{density(..., bw = "SJ-dpi")}.
#'
#' @return A list with components:
#' \describe{
#'   \item{bw}{Numeric vector of length 3: \code{c(epsilon, epsilon, delta)}.}
#'   \item{time}{Numeric vector of length \code{nrow(X)} giving temporal edge masses
#'     \eqn{\int_{tmin}^{tmax} K_\delta(u-t_i)\,du}.}
#'   \item{space}{Numeric vector of length \code{nrow(X)} giving spatial edge masses computed on the
#'     bounding box of \code{s.region}.}
#' }
#'
#' @details
#' The spatial window is converted to an \code{\link[spatstat.geom]{owin}} object, and the spatial
#' bandwidth is estimated using \code{bw.diggle} on the corresponding \code{\link[spatstat.geom]{ppp}}
#' object. The temporal bandwidth is taken as the bandwidth selected by \code{density} with
#' \code{bw = "SJ-dpi"} over the interval \code{t.region}.
#'
#' The returned edge-correction factors are kernel masses inside the window. If you use them for
#' intensity estimation with edge correction, typical usage is to divide by these factors.
#'
#' @references
#' Baddeley A, Rubak E, Turner R (2015).
#' \emph{Spatial Point Patterns: Methodology and Applications with R}.
#' Chapman and Hall/CRC Press.
#'
#' Diggle, P.J. (1985).
#' A Kernel Method for Smoothing Point Process Data.
#' \emph{Journal of the Royal Statistical Society, Series C},
#' 34, 138--147.
#'
#' @seealso \code{\link[spatstat.explore]{bw.diggle}}, \code{\link[stats]{density}}
#'
#' @examples
#'   set.seed(123)
#'   X <- cbind(runif(100), runif(100), runif(100, 0, 10))
#'   s.region <- matrix(c(0,0, 1,0, 1,1, 0,1), ncol = 2, byrow = TRUE)
#'   t.region <- c(0, 10)
#'   n.grid <- c(25, 25, 20)
#'   res <- calc.bandwidths.and.edgecorr(X, s.region, t.region, n.grid)
#'   str(res)
#'
#' @importFrom stats density pnorm
#' @importFrom spatstat.geom owin ppp boundingbox
#' @importFrom spatstat.explore bw.diggle
#' @export
calc.bandwidths.and.edgecorr <- function(X, s.region, t.region, n.grid,
                                         epsilon = NULL, delta = NULL) {
  # ---- basic checks (CRAN-friendly, informative) ----
  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) < 3L) {
    stop("`X` must be a numeric matrix/data.frame with at least 3 columns (x, y, t).")
  }
  if (anyNA(X[, 1:3, drop = FALSE])) stop("`X` contains NA in the first three columns.")
  if (!is.numeric(s.region) || !is.matrix(s.region) || ncol(s.region) != 2L || nrow(s.region) < 3L) {
    stop("`s.region` must be a numeric matrix with 2 columns and at least 3 rows (polygon vertices).")
  }
  if (!is.numeric(t.region) || length(t.region) != 2L || anyNA(t.region)) {
    stop("`t.region` must be a numeric vector of length 2: c(tmin, tmax).")
  }
  tmin <- min(t.region); tmax <- max(t.region)
  if (!(tmin < tmax)) stop("`t.region` must satisfy tmin < tmax.")
  if (!is.numeric(n.grid) || length(n.grid) < 3L || is.na(n.grid[3L]) || n.grid[3L] <= 1) {
    stop("`n.grid` must be a numeric/integer vector with length >= 3 and n.grid[3] > 1.")
  }
  n_t <- as.integer(n.grid[3L])

  # ---- construct window + point pattern ----
  ObsW <- spatstat.geom::owin(poly = list(x = s.region[, 1], y = s.region[, 2]))
  SPP  <- spatstat.geom::ppp(x = X[, 1], y = X[, 2], window = ObsW)

  # ---- spatial bandwidth ----
  if (is.null(epsilon)) {
    epsilon <- spatstat.explore::bw.diggle(SPP)
  }
  epsilon <- as.numeric(epsilon)
  if (!is.finite(epsilon) || epsilon <= 0) stop("`epsilon` must be a finite positive number.")

  # ---- temporal bandwidth ----
  if (is.null(delta)) {
    # density() may fail in degenerate cases; provide a simple fallback
    dens <- try(stats::density(X[, 3], bw = "SJ-dpi", n = n_t, from = tmin, to = tmax), silent = TRUE)
    if (inherits(dens, "try-error") || !is.finite(dens$bw) || dens$bw <= 0) {
      dens <- stats::density(X[, 3], bw = "nrd0", n = n_t, from = tmin, to = tmax)
    }
    delta <- dens$bw
  }
  delta <- as.numeric(delta)
  if (!is.finite(delta) || delta <= 0) stop("`delta` must be a finite positive number.")

  h <- c(epsilon, epsilon, delta)

  # ---- temporal edge mass (exact for Gaussian kernel on [tmin, tmax]) ----
  tvals <- X[, 3]
  edge.time <- stats::pnorm(tmax, mean = tvals, sd = delta) -
    stats::pnorm(tmin, mean = tvals, sd = delta)

  # ---- spatial edge mass using bounding-box approximation ----
  bb <- spatstat.geom::boundingbox(ObsW)
  xmin <- bb$xrange[1L]; xmax <- bb$xrange[2L]
  ymin <- bb$yrange[1L]; ymax <- bb$yrange[2L]

  xvals <- X[, 1]
  yvals <- X[, 2]

  xmass <- stats::pnorm(xmax, mean = xvals, sd = epsilon) -
    stats::pnorm(xmin, mean = xvals, sd = epsilon)
  ymass <- stats::pnorm(ymax, mean = yvals, sd = epsilon) -
    stats::pnorm(ymin, mean = yvals, sd = epsilon)

  edge.space <- xmass * ymass

  # Numerical safety: in extreme cases these can be ~0
  if (any(edge.time <= 0) || any(edge.space <= 0)) {
    warning("Some edge-correction masses are non-positive (numerical underflow or points far outside window).")
  }

  list(bw = h, time = edge.time, space = edge.space)
}
