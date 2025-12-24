#' Chi-squared test for first-order separability of a spatio-temporal point process
#'
#' Performs the classical (asymptotic) chi-squared test of first-order separability by
#' constructing a space--time contingency table of counts and applying a chi-squared test
#' of independence.
#'
#' The spatial domain is partitioned into \code{n.space} bins in each coordinate direction
#' (yielding \code{n.space^2} spatial cells), and the temporal domain is partitioned into
#' \code{n.time} intervals. Bin boundaries are defined using empirical quantiles of the
#' observed coordinates, with the first/last boundaries fixed to the provided spatial and
#' temporal windows.
#'
#' @param X A numeric matrix or data frame with at least three columns giving event coordinates
#'   \eqn{(x, y, t)}.
#' @param n.space Integer (>= 2). Number of bins per spatial axis. The contingency table has
#'   \code{n.space^2} rows.
#' @param n.time Integer (>= 2). Number of temporal bins (columns of the contingency table).
#' @param s.region Numeric vector of length 4 giving the spatial bounding box
#'   \code{c(xmin, xmax, ymin, ymax)}. Defaults to the unit box \code{c(0,1,0,1)}.
#' @param t.region Numeric vector of length 2 giving the temporal window \code{c(tmin, tmax)}
#'   with \code{tmin < tmax}. Defaults to \code{c(0,1)}.
#'
#' @return A list with components:
#' \describe{
#'   \item{chisq_s}{Numeric scalar. The chi-squared test statistic.}
#'   \item{chisq_p}{Numeric scalar. The p-value of the chi-squared test.}
#'   \item{counts}{Integer matrix of dimension \code{n.space^2} by \code{n.time} containing the space--time counts.}
#' }
#'
#' @details
#' This implementation uses \code{\link[stats]{chisq.test}} on the contingency table of space--time counts.
#' If expected counts are very small, the chi-squared approximation may be poor; in that case consider
#' using a Monte Carlo approach (e.g., block permutation) as implemented in \code{\link{chi2.test}}.
#'
#' @references
#' Ghorbani M., Vafaei N., Dvořák J., Myllymäki M. (2021).
#' Testing the first-order separability hypothesis for spatio-temporal point patterns.
#' \emph{Computational Statistics and Data Analysis}, \bold{161}, 107245.
#'
#' @seealso \code{\link{chi2.test}}, \code{\link[stats]{chisq.test}}
#'
#' @examples
#'
#' lambda <- get.lambda.function(N = 200, g = 50, model = 4)
#' Lmax <- get.lambda.max(N = 200, g = 50, model = 4)
#' X <- rstpoispp(lambda, Lmax)
#' result <- chisq.test.stPP(X, n.space = 2, n.time = 2)
#' print(result)
#'
#' @author
#' Jiří Dvořák \email{dvorak@karlin.mff.cuni.cz }\cr
#'
#' @importFrom stats chisq.test quantile
#' @export
chisq.test.stPP <- function(X,
                            n.space = 2L,
                            n.time  = 3L,
                            s.region = c(0, 1, 0, 1),
                            t.region = c(0, 1)) {

  # ---- validate inputs ----
  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) < 3L) stop("`X` must be numeric with at least 3 columns (x, y, t).")
  if (anyNA(X[, 1:3, drop = FALSE]) || any(!is.finite(X[, 1:3, drop = FALSE]))) {
    stop("The first three columns of `X` must be finite (no NA/NaN/Inf).")
  }

  n.space <- as.integer(n.space)
  n.time  <- as.integer(n.time)
  if (n.space < 2L) stop("`n.space` must be >= 2.")
  if (n.time < 2L)  stop("`n.time` must be >= 2.")

  if (!is.numeric(s.region) || length(s.region) != 4L || anyNA(s.region) || any(!is.finite(s.region))) {
    stop("`s.region` must be a finite numeric vector of length 4: c(xmin, xmax, ymin, ymax).")
  }
  xmin <- s.region[1L]; xmax <- s.region[2L]
  ymin <- s.region[3L]; ymax <- s.region[4L]
  if (!(xmin < xmax && ymin < ymax)) stop("`s.region` must satisfy xmin < xmax and ymin < ymax.")

  if (!is.numeric(t.region) || length(t.region) != 2L || anyNA(t.region) || any(!is.finite(t.region))) {
    stop("`t.region` must be a finite numeric vector of length 2: c(tmin, tmax).")
  }
  tmin <- t.region[1L]; tmax <- t.region[2L]
  if (!(tmin < tmax)) stop("`t.region` must satisfy tmin < tmax.")

  # ---- bin boundaries (quantiles, with fixed outer bounds) ----
  # Use type=7 default; keep as stats::quantile
  bx <- stats::quantile(X[, 1], probs = seq(0, 1, length.out = n.space + 1L), names = FALSE, type = 7)
  by <- stats::quantile(X[, 2], probs = seq(0, 1, length.out = n.space + 1L), names = FALSE, type = 7)
  bt <- stats::quantile(X[, 3], probs = seq(0, 1, length.out = n.time  + 1L), names = FALSE, type = 7)

  bx[1L] <- xmin; bx[n.space + 1L] <- xmax
  by[1L] <- ymin; by[n.space + 1L] <- ymax
  bt[1L] <- tmin; bt[n.time  + 1L] <- tmax

  # Ensure strictly increasing breaks (ties can happen with repeated coords)
  # If ties exist, chisq table may be degenerate; provide a clear error.
  if (any(diff(bx) <= 0) || any(diff(by) <= 0) || any(diff(bt) <= 0)) {
    stop("Non-increasing bin boundaries detected (likely due to ties). Reduce n.space/n.time or jitter data.")
  }

  # ---- assign bins using cut() (include lowest; right-closed by default) ----
  ix <- cut(X[, 1], breaks = bx, include.lowest = TRUE, labels = FALSE)
  iy <- cut(X[, 2], breaks = by, include.lowest = TRUE, labels = FALSE)
  it <- cut(X[, 3], breaks = bt, include.lowest = TRUE, labels = FALSE)

  # any points outside the specified windows become NA
  keep <- !is.na(ix) & !is.na(iy) & !is.na(it)
  if (!all(keep)) {
    # You can choose warning vs stop. Warning is often nicer.
    warning("Some events fall outside s.region/t.region and were ignored in the contingency table.")
    ix <- ix[keep]; iy <- iy[keep]; it <- it[keep]
  }

  # ---- build contingency table: rows = spatial cells, cols = time bins ----
  spatial_cell <- (ix - 1L) * n.space + iy  # 1..n.space^2
  counts <- matrix(0L, nrow = n.space^2L, ncol = n.time)
  for (k in seq_along(it)) {
    counts[spatial_cell[k], it[k]] <- counts[spatial_cell[k], it[k]] + 1L
  }

  # ---- chisq test ----
  tst <- stats::chisq.test(counts, correct = FALSE)

  list(chisq_s = unname(tst$statistic),
       chisq_p = unname(tst$p.value),
       counts  = counts)
}
