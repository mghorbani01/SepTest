#' Generate Permuted Versions of a Spatio-Temporal Point Pattern
#'
#' Implements two types of permutation procedures for resampling the time component of spatio-temporal point process data:
#' \describe{
#'   \item{\code{"pure"}}{Pure random permutation of the time coordinates.}
#'   \item{\code{"block"}}{Block permutation where the time dimension is divided into consecutive blocks, and permutations are applied at the block level.}
#' }
#' These procedures are used for generating surrogate datasets under the null hypothesis of first-order separability.
#'
#' @param X A numeric matrix or data frame with at least three columns, where the third column represents time.
#' @param nperm Integer. The number of permuted datasets to generate.
#' @param nblocks Integer. The number of temporal blocks to use for block permutation. Must be > 2.
#' @param method Character. The permutation strategy to use. One of \code{"pure"} or \code{"block"}.
#'
#' @return A list of \code{nperm} matrices. Each matrix is a permuted version of the original input \code{X}, where the third column (time) has been resampled based on the selected method.
#'
#' @importFrom combinat permn
#' @export
#'
#' @examples
#'
#' \donttest{
#' set.seed(123)
#' X <- cbind(runif(100), runif(100), sort(runif(100)))
#'
#' # Pure permutation
#' sims_pure <- sim.procedures(X, nperm = 10, method = "pure")
#' head(sims_pure[[1]])
#' # Block permutation
#' sims_block <- sim.procedures(X, nperm = 10, nblocks = 5, method = "block")
#'
#' # Visualize the first result from block permutation
#' plot_stpp(sims_block[[1]], type = "3D")
#' }

## 1) block and permutation
sim.procedures <- function(X, nperm = 1999, nblocks = 4, method = c("block", "pure")) {
  # Match method argument
  method <- match.arg(method)

  # Input checks
  if (!(is.matrix(X) || is.data.frame(X))) {
    stop("Input 'X' must be a matrix or a data frame.")
  }

  if (ncol(X) < 3) {
    stop("Input 'X' must have at least 3 columns (e.g., the third column is used for time or ordering).")
  }

  if (!is.numeric(nperm) || length(nperm) != 1 || nperm <= 0) {
    stop("Parameter 'nperm' must be a single positive numeric value.")
  }

  if (!is.numeric(nblocks) || length(nblocks) != 1 || nblocks <= 2) {
    stop("Parameter 'nblocks' (number of blocks) must be a single numeric value greater than 2.")
  }


  # Sort by time column (3rd column)
  X <- X[order(X[, 3]), ]

  # Handle "pure" permutation
  if (method == "pure") {
    sim <- replicate(nperm, {
       cbind(X[,1],X[,2], sample(X[, 3], replace = FALSE))
    }, simplify = FALSE)

    return(sim)
  }

  # Handle "block" permutation
  n <- trunc(nrow(X) / nblocks)
  size <- n * nblocks
  X_reduced <- X[1:size, ]

  # Create groups for block permutation
  groups <- cut(seq_len(size), breaks = nblocks, labels=FALSE)
  timeblocks <- split(X_reduced[, 3], groups)

  # Generate all or sampled permutations
  nall <- factorial(nblocks)
  t.permuted <- combinat::permn(1:nblocks)[-1] # Remove identity permutation

  if (nall > nperm) {
    t.permuted <- sample(t.permuted, nperm, replace=FALSE)
  }

  method <- match.arg(method)

  # Preallocate a list for permuted matrices
  len <- length(t.permuted)
  Y_mix <- vector("list", length=len)

  # Preallocate matrix outside the loop for faster updates
  for (i in seq_len(len)) {
    # Combine the permuted times into a new third column
    permuted_times <- do.call(c, timeblocks[t.permuted[[i]]])

    # Create a permuted matrix by updating the third column directly
    permuted_matrix <- X_reduced
    permuted_matrix[, 3] <- permuted_times

    # Add any remaining rows that were not permuted
    if (size < nrow(X)) {
      permuted_matrix <- rbind(permuted_matrix, X[(size + 1):nrow(X), ])
    }

    Y_mix[[i]] <- permuted_matrix
  }
  sim=Y_mix
  return(sim)
}

