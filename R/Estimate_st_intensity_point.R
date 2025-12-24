#' Kernel intensity estimates at observed spatio-temporal points
#'
#' Computes kernel-based spatial, temporal, separable, and non-separable intensity
#' estimates evaluated at the observed spatio-temporal event locations. The function
#' also returns the separability diagnostic \eqn{S_i} and global deviation measures
#' quantifying departures from first-order separability.
#'
#' @param X Numeric matrix/data.frame with three columns \code{(x,y,t)} giving event coordinates.
#' @param n.grid Integer. Included for API compatibility with grid-based routines; not used.
#' @param edge List with components \code{bw} (length 3), \code{space}, and \code{time}.
#'   \code{space} and \code{time} are Gaussian edge-correction masses evaluated at each event;
#'   each may be a scalar or a numeric vector of length \code{nrow(X)}.
#'
#' @details
#' Pairwise Gaussian kernel weights are computed in each dimension and diagonal
#' entries are set to zero to remove self-contributions.
#'
#' @return A list with components \code{S.fun}, deviation measures, and estimated
#' intensity components at the observed points.
#'
#' @seealso \code{\link[stats]{dnorm}}
#'
#' @examples
#' \donttest{
#' X <- cbind(stats::runif(50), stats::runif(50), stats::runif(50))
#' edge <- list(bw = c(0.1, 0.1, 0.1), space = 1, time = 1)
#' res <- estimate.intensity.point(X, n.grid = 50, edge = edge)
#' str(res)
#' }
#'
#' @importFrom stats dnorm
#' @export
estimate.intensity.point <- function(X, n.grid, edge) {

  # n.grid unused (kept for API compatibility)
  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) < 3L) stop("`X` must be a numeric matrix/data.frame with at least 3 columns.")
  X <- X[, 1:3, drop = FALSE]
  np <- nrow(X)
  if (np < 2L) stop("`X` must contain at least two points for pairwise estimation.")

  if (!is.list(edge) || is.null(edge$bw) || is.null(edge$space) || is.null(edge$time)) {
    stop("`edge` must be a list with components `bw`, `space`, and `time`.")
  }
  h <- edge$bw
  if (!is.numeric(h) || length(h) != 3L || any(!is.finite(h)) || any(h <= 0)) {
    stop("`edge$bw` must be a numeric vector of length 3 with strictly positive values.")
  }

  # allow scalar or length-n vectors; recycle scalars
  space_w <- edge$space
  time_w  <- edge$time

  if (!is.numeric(space_w) || any(!is.finite(space_w)) || any(space_w <= 0)) {
    stop("`edge$space` must be positive and finite (scalar or length nrow(X)).")
  }
  if (!is.numeric(time_w) || any(!is.finite(time_w)) || any(time_w <= 0)) {
    stop("`edge$time` must be positive and finite (scalar or length nrow(X)).")
  }

  if (length(space_w) == 1L) space_w <- rep(space_w, np)
  if (length(time_w)  == 1L) time_w  <- rep(time_w,  np)
  if (length(space_w) != np) stop("`edge$space` must have length 1 or nrow(X).")
  if (length(time_w)  != np) stop("`edge$time` must have length 1 or nrow(X).")

  # pairwise kernel contributions
  # divide row i by space_w[i] / time_w[i] (edge-correction masses at the data points)
  ax <- outer(X[, 1], X[, 1], stats::dnorm, sd = h[1L]) /
    matrix(rep(space_w, each = np), nrow = np, ncol = np, byrow = FALSE)

  ay <- outer(X[, 2], X[, 2], stats::dnorm, sd = h[2L])

  at <- outer(X[, 3], X[, 3], stats::dnorm, sd = h[3L]) /
    matrix(rep(time_w, each = np), nrow = np, ncol = np, byrow = FALSE)

  diag(ax) <- 0
  diag(ay) <- 0
  diag(at) <- 0

  Spat.intens <- rowSums(ax * t(ay))
  TeM.intens  <- rowSums(at)

  tiny <- .Machine$double.xmin
  Spat.intens <- pmax(Spat.intens, tiny)
  TeM.intens  <- pmax(TeM.intens,  tiny)

  sep.intens <- (Spat.intens * TeM.intens) / np
  sep.intens <- pmax(sep.intens, tiny)

  nonsep.intens <- rowSums(ax * t(ay) * at)
  nonsep.intens <- pmax(nonsep.intens, tiny)

  S.fun <- nonsep.intens / sep.intens

  deviation.t1 <- sum(abs(nonsep.intens - sep.intens))
  deviation.t2 <- sum(log(nonsep.intens) - log(sep.intens))
  deviation.t3 <- sum(S.fun)

  list(
    S.fun = S.fun,
    S.space = NULL,
    S.time = NULL,
    deviation.t1 = deviation.t1,
    deviation.t2 = deviation.t2,
    deviation.t3 = deviation.t3,
    SPat.intens = Spat.intens,
    TeM.intens = TeM.intens,
    nonsep.intens = nonsep.intens,
    sep.intens = sep.intens
  )
}
