#' Global envelope test for spatio-temporal separability using S-function
#'
#' Performs a global envelope test of the null hypothesis of first-order separability
#' for a spatio-temporal point process. The observed separability diagnostics
#' \eqn{S(u,t)}, \eqn{S_{\mathrm{space}}(u)} and/or \eqn{S_{\mathrm{time}}(t)} are compared
#' to a reference distribution obtained from permuted versions of the data.
#'
#' Two permutation strategies are supported:
#' \describe{
#'   \item{\code{"pure_per"}}{Pure permutation: randomly permutes the time coordinates.}
#'   \item{\code{"block_per"}}{Block permutation: permutes time in blocks to preserve short-range temporal dependence.}
#' }
#'
#' The \pkg{GET} package is used to construct global envelopes and compute p-values.
#'
#' @param X A numeric matrix or data frame with at least three columns giving \eqn{(x,y,t)}.
#' @param sim.procedure Character string specifying the permutation strategy:
#'   \code{"pure_per"} or \code{"block_per"}.
#' @param nsim Integer. Number of permutations for \code{"pure_per"}.
#' @param nblocks Integer (>= 2). Number of temporal blocks for block permutation.
#'   Used only for \code{"block_per"}.
#' @param nperm Integer. Number of block permutations for \code{"block_per"}.
#' @param n.grid Integer/numeric vector of length 3 specifying grid resolution in \eqn{x,y,t}.
#' @param s.region Numeric matrix with two columns specifying polygon vertices of the spatial window.
#' @param t.region Numeric vector of length 2 specifying temporal window \code{c(tmin,tmax)}.
#' @param owin Optional window of class \code{"owin"} (from \pkg{spatstat.geom}). If supplied,
#'   values outside the window are set to \code{NA} (pixels mode).
#' @param eps Optional numeric scalar (>0). Spatial bandwidth. If \code{NULL}, estimated internally.
#' @param del Optional numeric scalar (>0). Temporal bandwidth. If \code{NULL}, estimated internally.
#' @param tests Character vector indicating which diagnostics to test. Any of
#'   \code{"S.test"}, \code{"S.space.test"}, \code{"S.time.test"}.
#' @param GET.args Optional named list of extra arguments passed to
#'   \code{\link[GET]{global_envelope_test}} (e.g. \code{alternative}, \code{savefuns}, \code{nstep}).
#'
#' @details
#' The null hypothesis is
#' \deqn{H_0:\ \rho(u,t)=\rho_{\mathrm{space}}(u)\rho_{\mathrm{time}}(t).}
#' The function computes the chosen diagnostics using \code{\link{S.based.functions}}
#' on a pixel grid (\code{at="pixels"}) and applies \code{\link[GET]{global_envelope_test}}.
#'
#' To keep curve lengths identical (required by \pkg{GET}), any \code{NA} values induced by
#' an \code{owin} mask are removed using the **same indices** for the observed and all simulated curves.
#'
#' @return A list with components:
#' \describe{
#'   \item{Bandwidth_s}{Spatial bandwidth used.}
#'   \item{Bandwidth_t}{Temporal bandwidth used.}
#'   \item{S.test, S.space.test, S.time.test}{For each requested test: p-values for ERL and AREA envelopes,
#'   plus optional plots if created.}
#' }
#'
#' @author
#' Nafiseh Vafaei \email{nafiseh.vafaei@slu.se}\cr
#' Mohammad Ghorbani \email{mohammad.ghorbani@slu.se}
#'
#' @references
#' Ghorbani, M., Vafaei, N., Dvořák, J., and Myllymäki, M. (2021).
#' Testing the first-order separability hypothesis for spatio-temporal point patterns.
#' \emph{Computational Statistics & Data Analysis}, \bold{161}, 107245.
#'
#' @importFrom GET create_curve_set global_envelope_test
#' @seealso
#' \code{\link{S.based.functions}},
#' \code{\link{sim.procedures}},
#' \code{\link{block.permut}},
#' \code{\link[GET]{global_envelope_test}},
#' \code{\link[GET]{plot.global_envelope}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("GET", quietly = TRUE)) {
#'   set.seed(123)
#'   X <- cbind(stats::runif(100), stats::runif(100), stats::runif(100, 0, 1))
#'   s.region <- matrix(c(0,0, 1,0, 1,1, 0,1), ncol = 2, byrow = TRUE)
#'   t.region <- c(0, 1)
#'
#'   res <- global.envelope.test(
#'     X = X,
#'     sim.procedure = "pure_per",
#'     nsim = 19,
#'     n.grid = c(10,10,10),
#'     s.region = s.region,
#'     t.region = t.region,
#'     tests = c("S.test","S.time.test")
#'   )
#'   str(res)
#' }
#' }
#'
#' @export
global.envelope.test <- function(X,
                                 sim.procedure = c("pure_per", "block_per"),
                                 nsim = 25L,
                                 nblocks = 5L,
                                 nperm = 199L,
                                 n.grid = c(20L, 20L, 10L),
                                 s.region = matrix(c(0,0, 1,0, 1,1, 0,1), ncol = 2, byrow = TRUE),
                                 t.region = c(0, 1),
                                 owin = NULL,
                                 eps = NULL,
                                 del = NULL,
                                 tests = c("S.test", "S.space.test", "S.time.test"),
                                 GET.args = NULL) {

  if (!requireNamespace("GET", quietly = TRUE)) {
    stop("Package 'GET' is required for global.envelope.test().", call. = FALSE)
  }

  sim.procedure <- match.arg(sim.procedure)
  tests <- match.arg(tests, choices = c("S.test", "S.space.test", "S.time.test"), several.ok = TRUE)

  nsim <- as.integer(nsim)
  nblocks <- as.integer(nblocks)
  nperm <- as.integer(nperm)
  n.grid <- as.integer(n.grid)

  if (nsim < 1L) stop("`nsim` must be >= 1.", call. = FALSE)
  if (nblocks < 2L) stop("`nblocks` must be >= 2.", call. = FALSE)
  if (nperm < 1L) stop("`nperm` must be >= 1.", call. = FALSE)
  if (!is.null(GET.args)) {
    if (!is.list(GET.args) || is.null(names(GET.args)) || any(names(GET.args) == "")) {
      stop("`GET.args` must be a *named* list, or NULL.", call. = FALSE)
    }
  } else {
    GET.args <- list()
  }

  check.args(X, s.region, t.region, n.grid)

  # Validate owin only if supplied
  if (!is.null(owin)) {
    if (!requireNamespace("spatstat.geom", quietly = TRUE)) {
      stop("`owin` was supplied but package 'spatstat.geom' is not available.", call. = FALSE)
    }
    if (!inherits(owin, "owin")) stop("`owin` must be of class 'owin' (spatstat.geom).", call. = FALSE)
  }

  X <- as.matrix(X)[, 1:3, drop = FALSE]

  # --- observed diagnostics (pixels mode) ---
  obs_all <- S.based.functions(
    X = X, s.region = s.region, t.region = t.region,
    owin = owin, at = "pixels", n.grid = n.grid,
    epsilon = eps, delta = del, output = "all"
  )

  bandwidth_s <- obs_all$epsilon
  bandwidth_t <- obs_all$delta

  # helper: compute required objects for a given dataset, but only those requested
  compute_diag <- function(Xsim) {
    out <- list()
    if ("S.test" %in% tests)       out$S.fun   <- S.based.functions(Xsim, s.region, t.region, owin, at="pixels", n.grid=n.grid, epsilon=eps, delta=del, output="S.fun")$S
    if ("S.space.test" %in% tests) out$S.space <- S.based.functions(Xsim, s.region, t.region, owin, at="pixels", n.grid=n.grid, epsilon=eps, delta=del, output="S.space")$S
    if ("S.time.test" %in% tests)  out$S.time  <- S.based.functions(Xsim, s.region, t.region, owin, at="pixels", n.grid=n.grid, epsilon=eps, delta=del, output="S.time")$S
    out
  }

  # --- simulate datasets ---
  sim_list <- if (sim.procedure == "pure_per") {
    replicate(nsim, cbind(X[, 1:2], sample(X[, 3], replace = FALSE)), simplify = FALSE)
  } else {
    sim.procedures(X, nperm = nperm, nblocks = nblocks, method = "block")
  }

  sim_data <- lapply(sim_list, compute_diag)

  output_list <- list(Bandwidth_s = bandwidth_s, Bandwidth_t = bandwidth_t)

  # --- function to run GET safely with consistent NA removal ---
  run_get_test <- function(obs_curve, sim_curves, type) {
    obs_vec <- as.vector(obs_curve)

    # enforce identical indices removed across obs + all sims
    keep <- is.finite(obs_vec)
    if (!all(keep)) {
      obs_vec <- obs_vec[keep]
      sim_curves <- lapply(sim_curves, function(v) as.vector(v)[keep])
    } else {
      sim_curves <- lapply(sim_curves, function(v) as.vector(v))
    }

    sim_m <- do.call(cbind, sim_curves)  # length(obs) x nsim
    cset <- GET::create_curve_set(list(
      r = seq_along(obs_vec),
      obs = obs_vec,
      sim_m = sim_m
    ))

    gt <- do.call(GET::global_envelope_test, c(list(curve_set = cset, type = type), GET.args))
    list(
      p.value = attr(gt, "p"),
      envelope = gt
    )
  }

  # ---- S.test (S.fun) ----
  if ("S.test" %in% tests) {
    obs_curve <- obs_all$S.fun
    sim_curves <- lapply(sim_data, `[[`, "S.fun")

    erl  <- run_get_test(obs_curve, sim_curves, type = "erl")
    area <- run_get_test(obs_curve, sim_curves, type = "area")

    output_list$S.test <- list(
      p.value.erl  = erl$p.value,
      p.value.area = area$p.value,
      envelope.erl = erl$envelope,
      envelope.area = area$envelope
    )
  }

  # ---- S.space.test ----
  if ("S.space.test" %in% tests) {
    obs_curve <- obs_all$S.space
    sim_curves <- lapply(sim_data, `[[`, "S.space")

    erl  <- run_get_test(obs_curve, sim_curves, type = "erl")
    area <- run_get_test(obs_curve, sim_curves, type = "area")

    output_list$S.space.test <- list(
      p.value.erl  = erl$p.value,
      p.value.area = area$p.value,
      envelope.erl = erl$envelope,
      envelope.area = area$envelope
    )
  }

  # ---- S.time.test ----
  if ("S.time.test" %in% tests) {
    obs_curve <- obs_all$S.time
    sim_curves <- lapply(sim_data, `[[`, "S.time")

    erl  <- run_get_test(obs_curve, sim_curves, type = "erl")
    area <- run_get_test(obs_curve, sim_curves, type = "area")

    output_list$S.time.test <- list(
      p.value.erl  = erl$p.value,
      p.value.area = area$p.value,
      envelope.erl = erl$envelope,
      envelope.area = area$envelope
    )
  }

  output_list
}
