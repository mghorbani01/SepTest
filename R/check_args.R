#' Validate common arguments for spatio-temporal grid-based routines
#'
#' Checks the validity of common inputs used by spatio-temporal estimation and simulation
#' routines: event locations \code{X}, spatial window specification \code{s.region}, temporal window
#' \code{t.region}, and grid resolution \code{n.grid}.
#'
#' @param X A matrix or data frame with at least three columns giving event coordinates
#'   \eqn{(x, y, t)}. Only the first three columns are used for validation.
#' @param s.region Spatial observation region specification. Either:
#'   \itemize{
#'     \item a numeric vector of length 4 giving \code{c(xmin, xmax, ymin, ymax)}, or
#'     \item a numeric matrix with two columns giving polygon vertices \code{(x, y)} (at least 3 rows).
#'   }
#'   All values must be finite.
#' @param t.region Numeric vector of length 2 giving \code{c(tmin, tmax)} with \code{tmin < tmax}.
#'   Values must be finite.
#' @param n.grid Integer vector of length 3 giving the number of grid cells in the \eqn{x}, \eqn{y},
#'   and time directions. Must be positive.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of throwing an error if inputs are invalid.
#'
#' @details
#' This function is intended as a lightweight argument checker used internally by multiple functions.
#' It can also be called directly for debugging input issues.
#'
#' @examples
#'
#' X <- cbind(runif(100), runif(100), runif(100, 0, 10))
#' s.region <- matrix(c(0,0, 1,0, 1,1, 0,1), ncol = 2, byrow = TRUE)
#' t.region <- c(0, 10)
#' check.args(X, s.region, t.region, n.grid = c(25, 25, 20))
#'
#' @export

check.args <- function(X, s.region, t.region, n.grid = c(25L, 25L, 20L)) {

  # ---- X ----
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("`X` must be a matrix or data.frame.")
  }
  X <- as.matrix(X)
  if (ncol(X) < 3L) stop("`X` must have at least 3 columns (x, y, t).")
  if (!is.numeric(X[, 1:3, drop = FALSE])) stop("The first three columns of `X` must be numeric.")
  if (anyNA(X[, 1:3, drop = FALSE]) || any(!is.finite(X[, 1:3, drop = FALSE]))) {
    stop("All values in the first three columns of `X` must be finite (no NA/NaN/Inf).")
  }

  # ---- t.region ----
  if (!is.numeric(t.region) || length(t.region) != 2L) {
    stop("`t.region` must be a numeric vector of length 2: c(tmin, tmax).")
  }
  if (anyNA(t.region) || any(!is.finite(t.region))) stop("All values in `t.region` must be finite.")
  if (!(t.region[1L] < t.region[2L])) stop("`t.region` must satisfy tmin < tmax.")

  # ---- n.grid ----
  if (!is.numeric(n.grid) || length(n.grid) != 3L) {
    stop("`n.grid` must be a numeric/integer vector of length 3.")
  }
  if (anyNA(n.grid) || any(!is.finite(n.grid))) stop("All values in `n.grid` must be finite.")
  if (any(n.grid <= 0)) stop("All values in `n.grid` must be positive.")
  # grid cells should be integer-like
  if (any(abs(n.grid - round(n.grid)) > 0)) stop("`n.grid` must contain integer values.")
  n.grid <- as.integer(n.grid)

  # ---- s.region ----
  if (is.matrix(s.region)) {
    if (!is.numeric(s.region) || ncol(s.region) != 2L || nrow(s.region) < 3L) {
      stop("If `s.region` is a matrix, it must be numeric with 2 columns and at least 3 rows (polygon vertices).")
    }
    if (anyNA(s.region) || any(!is.finite(s.region))) stop("All values in `s.region` must be finite.")
  } else {
    if (!is.numeric(s.region) || length(s.region) != 4L) {
      stop("`s.region` must be either a 2-column vertex matrix or a numeric vector c(xmin, xmax, ymin, ymax).")
    }
    if (anyNA(s.region) || any(!is.finite(s.region))) stop("All values in `s.region` must be finite.")
    if (!(s.region[1L] < s.region[2L] && s.region[3L] < s.region[4L])) {
      stop("If `s.region` is a vector, it must satisfy xmin < xmax and ymin < ymax.")
    }
  }

  invisible(NULL)
}
