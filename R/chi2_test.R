#' Chi-squared test for first-order separability of a spatio-temporal point process
#'
#' Performs a chi-squared test for testing first-order separability of a spatio-temporal point process.
#' Two procedures are available:
#' \describe{
#'   \item{\code{"pure_per"}}{Classical asymptotic chi-squared test of independence on a space--time count table.}
#'   \item{\code{"block_per"}}{Monte Carlo permutation test based on block-wise permutations of the time component.}
#' }
#'
#' @param X A numeric matrix or data frame with at least three columns giving event coordinates
#'   \eqn{(x, y, t)}.
#' @param sim.procedure Character string specifying the procedure: \code{"pure_per"} or \code{"block_per"}.
#' @param nblocks Integer (>= 2). Number of temporal blocks used for block permutation (only for \code{"block_per"}).
#' @param nperm Integer (>= 1). Number of Monte Carlo permutations (only for \code{"block_per"}).
#' @param n.time Integer (>= 2). Number of temporal intervals in the contingency table.
#' @param n.space Integer (>= 2). The spatial domain is partitioned into \code{n.space} bins per axis
#'   (yielding \code{n.space^2} spatial cells) for the contingency table.
#' @param t.region Numeric vector of length 2 giving the temporal window \code{c(tmin, tmax)} with \code{tmin < tmax}.
#' @param s.region Spatial window specification. By default, the bounding box \code{c(0, 1, 0, 1)}
#'   corresponding to \code{c(xmin, xmax, ymin, ymax)}. Passed to \code{\link{chisq.test.stPP}}.
#'
#' @return Numeric scalar: the p-value of the test.
#'
#' @details
#' The classical procedure (\code{"pure_per"}) applies a chi-squared test of independence to the
#' \code{n.space^2} by \code{n.time} contingency table of counts.
#'
#' The permutation procedure (\code{"block_per"}) generates \code{nperm} block-permuted datasets under the null
#' using \code{\link{sim.procedures}} with \code{method = "block"}, recomputes the chi-squared statistic for each,
#' and returns a Monte Carlo p-value computed as \eqn{(1 + \#\{T_i \ge T_{obs}\})/(nperm + 1)}.
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
#' Ghorbani, M., Vafaei, N. and Myllymäki, M. (2025).
#' A kernel-based test for the first-order separability of spatio-temporal point processes,
#' \emph{TEST}.
#'
#' @seealso \code{\link{chisq.test.stPP}}, \code{\link{sim.procedures}}, \code{\link{block.permut}}
#'
#' @examples
#'
#' set.seed(124)
#' lambda <- get.lambda.function(N = 200, g = 50, model = 4)
#' Lmax <- get.lambda.max(N = 200, g = 50, model = 4)
#' X <- rstpoispp(lambda, Lmax)
#'
#'
#' # Classical chi-squared test
#' chi2.test(X, sim.procedure = "pure_per", n.time = 2, n.space = 3)
#'
#' # Monte Carlo permutation test with blocks
#' chi2.test(X, sim.procedure = "block_per", nblocks = 5, nperm = 100)
#'
#' @export

chi2.test <- function(X,
                      sim.procedure = c("pure_per", "block_per"),
                      nblocks = 5L, nperm = 199L,
                      n.time = 2L, n.space = 3L,
                      t.region = c(0, 1),
                      s.region = c(0, 1, 0, 1)) {

  sim.procedure <- match.arg(sim.procedure)

  # ---- validate inputs ----
  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) < 3L) stop("`X` must be numeric with at least 3 columns (x, y, t).")
  if (anyNA(X[, 1:3, drop = FALSE]) || any(!is.finite(X[, 1:3, drop = FALSE]))) {
    stop("The first three columns of `X` must be finite (no NA/NaN/Inf).")
  }

  if (!is.numeric(t.region) || length(t.region) != 2L || anyNA(t.region) || any(!is.finite(t.region))) {
    stop("`t.region` must be a finite numeric vector of length 2.")
  }
  if (!(t.region[1L] < t.region[2L])) stop("`t.region` must satisfy tmin < tmax.")

  # s.region: expect bbox c(xmin,xmax,ymin,ymax) (simple and checkable)
  if (!is.numeric(s.region) || length(s.region) != 4L || anyNA(s.region) || any(!is.finite(s.region))) {
    stop("`s.region` must be a finite numeric vector of length 4: c(xmin, xmax, ymin, ymax).")
  }
  if (!(s.region[1L] < s.region[2L] && s.region[3L] < s.region[4L])) {
    stop("`s.region` must satisfy xmin < xmax and ymin < ymax.")
  }

  n.time  <- as.integer(n.time)
  n.space <- as.integer(n.space)
  if (n.time < 2L) stop("`n.time` must be >= 2.")
  if (n.space < 2L) stop("`n.space` must be >= 2.")

  if (sim.procedure == "pure_per") {
    res <- chisq.test.stPP(X = X, n.space = n.space, n.time = n.time,
                           s.region = s.region, t.region = t.region)
    return(res$chisq_p)
  }

  # ---- block permutation procedure ----
  nblocks <- as.integer(nblocks)
  nperm   <- as.integer(nperm)
  if (nblocks < 2L) stop("`nblocks` must be >= 2 for block permutation.")
  if (nperm < 1L) stop("`nperm` must be >= 1 for block permutation.")

  sims <- sim.procedures(X, nperm = nperm, nblocks = nblocks, method = "block")

  obs <- chisq.test.stPP(X = X, n.space = n.space, n.time = n.time,
                         s.region = s.region, t.region = t.region)$chisq_s

  sim_stats <- vapply(sims, function(Xsim) {
    chisq.test.stPP(X = Xsim, n.space = n.space, n.time = n.time,
                    s.region = s.region, t.region = t.region)$chisq_s
  }, numeric(1))

  # Monte Carlo p-value with +1 correction
  (1 + sum(sim_stats >= obs)) / (length(sim_stats) + 1)
}
