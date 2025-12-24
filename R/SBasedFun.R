#' Compute S-based test function for testing the null hypothesis of first-order separability
#'
#' Computes kernel-based estimates of the spatio-temporal intensity and related
#' separability diagnostics, either on a regular spatio-temporal grid (\code{"pixels"})
#' or at the observed event locations (\code{"points"}).
#'
#' The function is a wrapper that (i) validates inputs, (ii) computes bandwidths and
#' Gaussian edge-correction masses via \code{\link{calc.bandwidths.and.edgecorr}}, and
#' (iii) delegates the actual estimation to \code{\link{estimate.intensity.pixel}} or
#' \code{\link{estimate.intensity.point}}.
#'
#' @param X Numeric matrix/data.frame with three columns giving \eqn{(x,y,t)}.
#' @param s.region Numeric matrix with two columns giving polygon vertices of the spatial window.
#' @param t.region Numeric vector of length 2 giving the temporal window \code{c(tmin, tmax)}.
#' @param owin Optional spatial window of class \code{"owin"} (from \pkg{spatstat.geom}).
#'   Used only when \code{at="pixels"} to set values outside the window to \code{NA}.
#' @param at Character string: \code{"pixels"} or \code{"points"}.
#' @param n.grid Integer vector of length 3 giving the grid resolution in \eqn{x}, \eqn{y}, and \eqn{t}
#'   (used when \code{at="pixels"}; still checked for length/positivity in both modes).
#' @param epsilon Optional numeric scalar (>0). Spatial bandwidth. If \code{NULL}, estimated internally.
#' @param delta Optional numeric scalar (>0). Temporal bandwidth. If \code{NULL}, estimated internally.
#' @param output Character string selecting which component to return. Use \code{"all"} (default) to return
#'   all available components; otherwise return only the selected component along with bandwidths and coordinates.
#'
#' @details
#' When \code{at="points"}, spatial/temporal profiles such as \code{S.space} and \code{S.time} are typically
#' not defined and are returned as \code{NULL} by the pointwise routine.
#'
#' @return
#' If \code{output="all"}, returns the full list produced by the chosen computation routine, augmented with
#' \code{epsilon} and \code{delta}.
#'
#' If \code{output} is a single component name, returns a list with:
#' \describe{
#'   \item{S}{The requested component.}
#'   \item{epsilon}{Spatial bandwidth used.}
#'   \item{delta}{Temporal bandwidth used.}
#'   \item{x,y,t}{Coordinates returned by the underlying routine (may be \code{NULL}).}
#' }
#'
#' @examples
#' \donttest{
#' X <- cbind(stats::runif(50), stats::runif(50), stats::runif(50))
#' s.region <- matrix(c(0,0, 1,0, 1,1, 0,1), ncol = 2, byrow = TRUE)
#' t.region <- c(0, 1)
#' res_all <- S.based.functions(X, s.region, t.region, at = "points")
#' res_Sfun <- S.based.functions(X, s.region, t.region, at = "points", output = "S.fun")
#' }
#'
#' @export
S.based.functions <- function(X, s.region, t.region,
                              owin = NULL,
                              at = c("pixels", "points"),
                              n.grid = c(25L, 25L, 20L),
                              epsilon = NULL, delta = NULL,
                              output = "all") {

  at <- match.arg(at)

  # ---- validate basics (reuse your checker) ----
  check.args(X, s.region, t.region, n.grid)

  # coerce and validate bandwidth inputs (CRAN-clean)
  if (!is.null(epsilon)) {
    if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) {
      stop("`epsilon` must be a single finite positive numeric value, or NULL.")
    }
  }
  if (!is.null(delta)) {
    if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta <= 0) {
      stop("`delta` must be a single finite positive numeric value, or NULL.")
    }
  }
  if (!is.character(output) || length(output) != 1L) stop("`output` must be a single character string.")

  # owin validity (only matters for pixels)
  if (!is.null(owin) && at == "pixels") {
    if (!requireNamespace("spatstat.geom", quietly = TRUE)) {
      stop("`owin` was supplied but package 'spatstat.geom' is not available.", call. = FALSE)
    }
    if (!inherits(owin, "owin")) stop("`owin` must be of class 'owin' (spatstat.geom).")
  }


  # ---- bandwidths + edge correction masses ----
  # NOTE: consistent with your earlier definition: no `owin` argument here.
  edge <- calc.bandwidths.and.edgecorr(X = X, s.region = s.region, t.region = t.region,
                                       n.grid = n.grid, epsilon = epsilon, delta = delta)

  # ---- compute requested mode ----
  result <- switch(
    at,
    pixels = estimate.intensity.pixel(X, s.region, t.region, n.grid, edge),
    points = estimate.intensity.point(X, n.grid, edge)
  )

  # always attach bandwidths
  result$epsilon <- edge$bw[1L]
  result$delta <- edge$bw[3L]

  # output selection
  if (identical(output, "all")) {
    return(result)
  }

  if (!output %in% names(result)) {
    stop("Invalid `output`. Choose one of: ", paste(names(result), collapse = ", "))
  }

  if (at == "points" && output %in% c("S.space", "S.time")) {
    warning(sprintf("%s is typically not computed in 'points' mode; returning the stored value (may be NULL).", output))
  }

  list(
    S = result[[output]],
    epsilon = result$epsilon,
    delta = result$delta,
    x = result$x, y = result$y, t = result$t
  )
}
