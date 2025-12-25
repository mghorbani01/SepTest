#' Kernel-based intensity estimation on a space-time grid and its components,
#' and test statistics for first-order separability
#'
#' Computes kernel-smoothed estimates of spatial, temporal, separable, and non-separable
#' spatio-temporal intensity functions on a regular space-time grid, together with
#' separability diagnostics used in first-order separability analysis.
#'
#' The estimator uses product Gaussian kernels with supplied bandwidths and (Gaussian)
#' edge-correction factors, typically produced by \code{\link{calc.bandwidths.and.edgecorr}}.
#'
#' @param X Numeric matrix/data.frame with three columns \code{(x, y, t)} giving observed events.
#' @param s.region Numeric matrix with two columns defining the spatial window (typically polygon vertices).
#'   Grid limits are taken as \code{range(s.region[,1])} and \code{range(s.region[,2])}.
#' @param t.region Numeric vector of length 2 giving the temporal window \code{c(tmin, tmax)}.
#' @param n.grid Integer vector of length 3 giving grid resolution in \eqn{x}, \eqn{y}, and \eqn{t}.
#' @param edge List with components \code{bw} (length 3), \code{space} (length \code{nrow(X)}),
#'   and \code{time} (length \code{nrow(X)}).
#' @param owin Optional observation window of class \code{"owin"} (from \pkg{spatstat.geom}).
#'   If provided, intensity estimates outside the window are set to \code{NA}.
#'
#' @return A list with grid coordinates \code{x,y,t}, intensity estimates, the diagnostic \code{S.fun},
#'   its marginal summaries \code{S.space} and \code{S.time}, and deviation measures.
#'
#' @references
#' Ghorbani, M., Vafaei, N., Dvořák, J., and Myllymäki, M. (2021).
#' Testing the first-order separability hypothesis for spatio-temporal point patterns.
#' \emph{Computational Statistics & Data Analysis}, 161, 107245.
#'
#' @seealso \code{\link{S.based.functions}}, \code{\link{calc.bandwidths.and.edgecorr}}
#'
#' @examples
#' \donttest{
#' n <- 100
#' X <- cbind(x = stats::runif(n), y = stats::runif(n), t = stats::runif(n, 0, 10))
#' s.region <- matrix(c(0,0, 1,0, 1,1, 0,1), ncol=2, byrow=TRUE)
#' t.region <- c(0, 10)
#' n.grid <- c(10, 10, 5)
#' edge <- list(bw = c(0.05, 0.05, 0.5), space = rep(1, n), time = rep(1, n))
#' res <- estimate.intensity.pixel(X, s.region, t.region, n.grid, edge)
#' str(res)
#' }
#'
#' @importFrom stats dnorm
#' @export
estimate.intensity.pixel <- function(X, s.region, t.region, n.grid, edge, owin = NULL) {

  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) < 3L) stop("`X` must be a numeric matrix/data.frame with at least 3 columns.")
  X <- X[, 1:3, drop = FALSE]
  np <- nrow(X)
  if (np < 1L) stop("`X` must contain at least one event.")

  if (!is.numeric(s.region) || !is.matrix(s.region) || ncol(s.region) != 2L) {
    stop("`s.region` must be a numeric matrix with 2 columns.")
  }
  if (!is.numeric(t.region) || length(t.region) != 2L || !(t.region[1] < t.region[2])) {
    stop("`t.region` must be a numeric vector c(tmin,tmax) with tmin < tmax.")
  }
  if (!is.numeric(n.grid) || length(n.grid) != 3L || any(!is.finite(n.grid)) || any(n.grid < 2)) {
    stop("`n.grid` must be an integer/numeric vector of length 3 with all entries >= 2.")
  }
  n.grid <- as.integer(n.grid)

  if (!is.list(edge) || is.null(edge$bw) || is.null(edge$space) || is.null(edge$time)) {
    stop("`edge` must be a list with components `bw`, `space`, and `time`.")
  }
  h <- edge$bw
  if (!is.numeric(h) || length(h) != 3L || any(!is.finite(h)) || any(h <= 0)) {
    stop("`edge$bw` must be a numeric vector of length 3 with strictly positive values.")
  }
  if (!is.numeric(edge$space) || length(edge$space) != np || any(!is.finite(edge$space)) || any(edge$space <= 0)) {
    stop("`edge$space` must be a positive numeric vector of length nrow(X).")
  }
  if (!is.numeric(edge$time) || length(edge$time) != np || any(!is.finite(edge$time)) || any(edge$time <= 0)) {
    stop("`edge$time` must be a positive numeric vector of length nrow(X).")
  }

  use_mask <- !is.null(owin)
  if (use_mask) {
    if (!requireNamespace("spatstat.geom", quietly = TRUE)) {
      stop("`owin` was provided but package 'spatstat.geom' is not available.", call. = FALSE)
    }
    if (!inherits(owin, "owin")) stop("`owin` must be of class 'owin' (spatstat.geom).")
  }

  xmin <- min(s.region[, 1]); xmax <- max(s.region[, 1])
  ymin <- min(s.region[, 2]); ymax <- max(s.region[, 2])
  tmin <- t.region[1];        tmax <- t.region[2]

  dx <- (xmax - xmin) / n.grid[1L]
  dy <- (ymax - ymin) / n.grid[2L]
  dt <- (tmax - tmin) / n.grid[3L]

  gx <- seq.int(xmin + dx/2, xmax - dx/2, length.out = n.grid[1L])
  gy <- seq.int(ymin + dy/2, ymax - dy/2, length.out = n.grid[2L])
  gt <- seq.int(tmin + dt/2, tmax - dt/2, length.out = n.grid[3L])

  # contributions
  ax <- outer(gx, X[, 1], stats::dnorm, sd = h[1L]) /
    matrix(rep(edge$space, each = n.grid[1L]), nrow = n.grid[1L], ncol = np, byrow = FALSE)

  ay <- outer(gy, X[, 2], stats::dnorm, sd = h[2L])

  at <- outer(gt, X[, 3], stats::dnorm, sd = h[3L]) /
    matrix(rep(edge$time, each = n.grid[3L]), nrow = n.grid[3L], ncol = np, byrow = FALSE)

  # spatial intensity (x by y)
  Spat.intens <- tcrossprod(ax, ay)

  tiny <- .Machine$double.xmin
  Spat.intens <- pmax(Spat.intens, tiny)

  maskY <- NULL
  if (use_mask) {
    m <- spatstat.geom::as.mask(owin, xy = list(x = gx, y = gy))
    maskY <- m$m
    # ensure mask orientation matches Spat.intens
    if (!identical(dim(maskY), dim(Spat.intens))) maskY <- t(maskY)
    Spat.intens[!maskY] <- NA_real_
  }

  # temporal intensity (t)
  TeM.intens <- rowSums(at)

  # separable intensity (x by y by t)
  SEP.intensity <- (Spat.intens %o% TeM.intens) / np
  SEP.intensity <- pmax(SEP.intensity, tiny)

  # non-separable intensity
  non.sep <- array(0, dim = n.grid)
  t.ay <- t(ay)  # (np x ny)
  for (k in seq_len(n.grid[3L])) {
    t.ay.tk <- t.ay * at[k, ]
    non.sep[, , k] <- ax %*% t.ay.tk
    if (use_mask) non.sep[, , k][!maskY] <- NA_real_
  }
  non.sep <- pmax(non.sep, tiny)

  # diagnostics
  g.s <- dx * dy
  t.s <- dt

  S.fun <- non.sep / SEP.intensity
  S.space <- apply(S.fun, 1:2, sum, na.rm = TRUE) * t.s
  S.time  <- apply(S.fun, 3,   sum, na.rm = TRUE) * g.s

  deviation.t1 <- sum(abs(non.sep - SEP.intensity), na.rm = TRUE) * t.s * g.s
  deviation.t3 <- sum(log(non.sep) - log(SEP.intensity), na.rm = TRUE) * t.s * g.s
  deviation.t4 <- sum(S.fun, na.rm = TRUE) * t.s * g.s

  list(
    x = gx, y = gy, t = gt,
    S.fun = S.fun, S.space = S.space, S.time = S.time,
    deviation.t1 = deviation.t1, deviation.t3 = deviation.t3, deviation.t4 = deviation.t4,
    SPat.intens = Spat.intens, TeM.intens = TeM.intens,
    nonsep.intens = non.sep, sep.intens = SEP.intensity
  )
}
