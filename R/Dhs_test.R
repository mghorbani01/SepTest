#' dHSIC test for first-order separability of a spatio-temporal point process
#'
#' Performs a nonparametric test of first-order separability between space and time in a
#' spatio-temporal point process using the distance-based Hilbert--Schmidt independence
#' criterion (dHSIC). The test statistic evaluates whether the spatio-temporal intensity can be factorized as a function of
#' spatial coordinate times as a function of temporal coordinate.
#'
#' Two permutation strategies are supported:
#' \describe{
#'   \item{\code{"pure_per"}}{Randomly permutes the time component across events.}
#'   \item{\code{"block_per"}}{Uses block-wise permutation of the time component via \code{\link{sim.procedures}}
#'     with \code{method = "block"} to preserve short-range temporal dependence.}
#' }
#'
#' @param X A numeric matrix or data frame with at least three columns giving event coordinates
#'   \eqn{(x, y, t)}.
#' @param sim.procedure Character string specifying the permutation strategy:
#'   \code{"pure_per"} or \code{"block_per"}.
#' @param nblocks Integer (>= 2). Number of temporal blocks for block permutation (only for \code{"block_per"}).
#' @param nperm Integer (>= 1). Number of block permutations (only for \code{"block_per"}).
#' @param nsim Integer (>= 1). Number of pure permutations (only for \code{"pure_per"}).
#' @param bandwidth Optional numeric. Fixed bandwidth to use with \code{kernel = "gaussian.fixed"} in \code{dHSIC::dhsic}.
#'   If provided, the function also returns a fixed-bandwidth Monte Carlo p-value. If \code{NULL},
#'   only the adaptive-bandwidth Gaussian kernel (\code{kernel = "gaussian"}) is used.
#'
#' @return A list with components:
#' \describe{
#'   \item{p.value}{Monte Carlo p-value based on the adaptive-bandwidth Gaussian kernel.}
#'   \item{p.value.bw}{Monte Carlo p-value based on the fixed-bandwidth Gaussian kernel, or \code{NA} if \code{bandwidth = NULL}.}
#'   \item{bandwidth_data}{Bandwidth selected by \code{dHSIC::dhsic(..., kernel = "gaussian")}.}
#' }
#'
#' @details
#' The Monte Carlo p-value is computed with the standard +1 correction:
#' \eqn{(1 + \#\{T_i \ge T_{obs}\})/(B + 1)}, where \eqn{B} is the number of permutations.
#'
#' @references
#' Ghorbani, M., Vafaei, N. and Myllymäki, M. (2025).
#' A kernel-based test for the first-order separability of spatio-temporal point processes,
#' \emph{TEST}.
#'
#' @seealso \code{\link{sim.procedures}}, \code{\link{block.permut}}, \code{\link{chi2.test}}
#'
#' @examples
#' if (requireNamespace("dHSIC", quietly = TRUE)) {
#'   set.seed(123)
#'   X <- cbind(runif(100), runif(100), runif(100, 0, 10))
#'
#'   # Pure permutation test
#'   result <- dHS.test(sim.procedure = "pure_per",
#'                      X = X, nsim = 199, bandwidth = 0.05)
#'   print(result$p.value)
#'
#'   # Block permutation test
#'   result_block <- dHS.test(sim.procedure = "block_per", X = X,
#'                            nblocks = 5, nperm = 100, bandwidth = 0.05)
#'   print(result_block$p.value.bw)
#' }
#' @note
#' The dHSIC method is implemented via the \pkg{dHSIC} package. When \code{sim.procedure = "pure_per"},
#' \code{dHS.test()} internally calls \code{\link{global.envelope.test}} for computational efficiency.
#'
#' @importFrom dHSIC dhsic
#'
#' @export

dHS.test <- function(X,
                     sim.procedure = c("pure_per", "block_per"),
                     nblocks = 7L,
                     nperm = 1999L,
                     nsim = 199L,
                     bandwidth = NULL) {

  sim.procedure <- match.arg(sim.procedure)

  # ---- dependency check (CRAN-friendly) ----
  if (!requireNamespace("dHSIC", quietly = TRUE)) {
    stop("Package 'dHSIC' is required for dHS.test(). Please install it.", call. = FALSE)
  }

  # ---- validate inputs ----
  X <- as.matrix(X)
  if (!is.numeric(X) || ncol(X) < 3L) stop("`X` must be numeric with at least 3 columns (x, y, t).")
  if (anyNA(X[, 1:3, drop = FALSE]) || any(!is.finite(X[, 1:3, drop = FALSE]))) {
    stop("The first three columns of `X` must be finite (no NA/NaN/Inf).")
  }

  if (!is.null(bandwidth)) {
    if (!is.numeric(bandwidth) || length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0) {
      stop("`bandwidth` must be a single finite positive numeric value, or NULL.")
    }
    bandwidth <- as.numeric(bandwidth)
  }

  # observed statistic (adaptive gaussian)
  res_obs <- dHSIC::dhsic(X[, 1:2, drop = FALSE], X[, 3], kernel = "gaussian")
  obs_stat <- res_obs$dHSIC
  bw_hat <- res_obs$bandwidth

  # observed statistic (fixed bw), optional
  obs_stat_bw <- NA_real_
  if (!is.null(bandwidth)) {
    obs_stat_bw <- dHSIC::dhsic(X[, 1:2, drop = FALSE], X[, 3],
                                kernel = "gaussian.fixed", bandwidth = bandwidth)$dHSIC
  }

  # ---- generate null statistics ----
  if (sim.procedure == "pure_per") {

    nsim <- as.integer(nsim)
    if (nsim < 1L) stop("`nsim` must be >= 1 for sim.procedure = 'pure_per'.")

    sim_stats <- numeric(nsim)
    sim_stats_bw <- if (!is.null(bandwidth)) numeric(nsim) else NULL

    for (b in seq_len(nsim)) {
      Xb_t <- sample(X[, 3], replace = FALSE)
      sim_stats[b] <- dHSIC::dhsic(X[, 1:2, drop = FALSE], Xb_t, kernel = "gaussian")$dHSIC
      if (!is.null(bandwidth)) {
        sim_stats_bw[b] <- dHSIC::dhsic(X[, 1:2, drop = FALSE], Xb_t,
                                        kernel = "gaussian.fixed", bandwidth = bandwidth)$dHSIC
      }
    }

  } else { # block_per

    nblocks <- as.integer(nblocks)
    nperm   <- as.integer(nperm)
    if (nblocks < 2L) stop("`nblocks` must be >= 2 for sim.procedure = 'block_per'.")
    if (nperm < 1L)   stop("`nperm` must be >= 1 for sim.procedure = 'block_per'.")

    sims <- sim.procedures(X, nperm = nperm, nblocks = nblocks, method = "block")

    sim_stats <- vapply(sims, function(Xsim) {
      dHSIC::dhsic(Xsim[, 1:2, drop = FALSE], Xsim[, 3], kernel = "gaussian")$dHSIC
    }, numeric(1))

    sim_stats_bw <- NULL
    if (!is.null(bandwidth)) {
      sim_stats_bw <- vapply(sims, function(Xsim) {
        dHSIC::dhsic(Xsim[, 1:2, drop = FALSE], Xsim[, 3],
                     kernel = "gaussian.fixed", bandwidth = bandwidth)$dHSIC
      }, numeric(1))
    }
  }

  # ---- Monte Carlo p-values (+1 correction) ----
  p.value <- (1 + sum(sim_stats >= obs_stat)) / (length(sim_stats) + 1)

  p.value.bw <- NA_real_
  if (!is.null(bandwidth)) {
    p.value.bw <- (1 + sum(sim_stats_bw >= obs_stat_bw)) / (length(sim_stats_bw) + 1)
  }

  list(
    p.value = p.value,
    p.value.bw = p.value.bw,
    bandwidth_data = bw_hat
  )
}
