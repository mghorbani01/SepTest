#' Apply a circular random shift to the temporal component of a spatio-temporal point pattern
#'
#' Performs a circular random shift of the temporal coordinate in a spatio-temporal
#' point pattern. This operation preserves the spatial configuration while
#' randomizing the temporal component under the assumption of temporal stationarity.
#' The shift amount is drawn uniformly from \eqn{[0, 1]} and applied modulo 1, ensuring
#' that the time window is maintained. The resulting dataset can be used to construct
#' null models for hypothesis testing of first-order separability or temporal independence.
#'
#' @param X A numeric matrix or data frame containing the spatio-temporal point pattern.
#' Must include at least one numeric column representing time.
#' @param shifted_col Integer index specifying which column to shift (typically the
#' time coordinate). Default is \code{3}.
#'
#' @details
#' The circular random shift is a common resampling procedure for generating
#' null models of temporal randomness while preserving the overall temporal
#' marginal distribution and spatial structure.
#'
#' For each dataset, a single uniform random shift value
#' \eqn{\Delta \sim \mathrm{Uniform}(0,1)} is drawn and added to the temporal
#' coordinate. The shifted times are then wrapped around the unit interval:
#'
#' \deqn{
#' t_i^{*} = (t_i + \Delta) \bmod 1, \quad i = 1, \dots, n.
#' }
#'
#' @return A matrix or data frame (matching the input type) with the time column shifted modulo 1.
#' @importFrom stats runif
#'
#' @references
#' Ghorbani, M., Vafaei, N., Dvořák, J., and Myllymäki, M. (2021).
#' Testing the first-order separability hypothesis for spatio-temporal point patterns.
#' \emph{Computational Statistics & Data Analysis}, 161, 107245.
#' @export
#'
#' @examples
#'
#' \donttest{
#' set.seed(123)
#' X <- cbind(runif(100), runif(100), runif(100))  # x, y, t
#' X_shifted <- random.shift(X)
#'
#' # Compare original and shifted time values
#' head(X[,3])
#' head(X_shifted[,3])
## 2) Random shift
#'
#' # Verify shift visually
#' plot(X[, 3], type = "o", col = "blue", ylab = "Time", xlab = "Index",
#'      main = "Original vs Shifted Times")
#' lines(X_shifted[, 3], type = "o", col = "red")
#' legend("topright", legend = c("Original", "Shifted"),
#'        col = c("blue", "red"), lty = 1, pch = 1)
#'  }
random.shift <- function(X, shifted_col = 3) {
  # Input validation
  if (!(is.matrix(X) || is.data.frame(X))) {
    stop("Input 'X' must be a matrix or a data frame.")
  }
  if (!is.numeric(shifted_col) || shifted_col < 1 || shifted_col > ncol(X)) {
    stop("Invalid time column index.")
  }

  # Extract and shift the time column
  t <- X[, shifted_col]
  if (!is.numeric(t)) {
    stop("Time column must be numeric.")
  }

  shift <- runif(1, 0, 1)
  t_new <- (t + shift) %% 1

  # Return a copy of X with shifted time column
  X_new <- X
  X_new[, shifted_col] <- t_new

  return(X_new)
}

