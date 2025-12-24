#' Block permutation of the temporal component in a spatio-temporal point pattern
#'
#' Permutes the temporal component of a spatio-temporal dataset in a block-wise manner while keeping
#' the spatial coordinates fixed. This is used to generate permuted replicates under the null model
#' of first-order separability.
#'
#' @param nblocks Integer (>= 2). Number of consecutive temporal blocks after ordering events by time.
#' @param X Numeric matrix or data frame with at least three columns \code{(x, y, t)}. Each row displays one event.
#'   The third column is interpreted as the time coordinate.
#' @param nperm Integer (>= 1). Number of permuted datasets to generate. At most
#'   \code{factorial(nblocks) - 1} distinct non-identity block permutations exist.
#'
#' @details The function first orders the events by time and partitions the ordered sequence into
#' \code{nblocks} consecutive blocks of equal size. The block labels are permuted (excluding the identity
#' permutation), and the time values are reassigned according to the permuted block order.
#'
#' If \code{nrow(X)} is not divisible by \code{nblocks}, the last \code{nrow(X) %% nblocks} events are not
#' included in the block permutation and are appended unchanged to each permuted dataset.
#'
#' For details of the block permutation procedure, see the Supplementary Materials in
#' Ghorbani et al. (2025).
#'
#'\code{\link{sim.procedures}}  covers both pure and block permutation methods.
#'
#' @return A list of length \code{min(nperm, factorial(nblocks) - 1)}. Each element is a matrix with the same
#' number of columns as \code{X}; the third column contains the block-permuted time values.
#'
#' @seealso \code{\link{sim.procedures}}
#'
#' @references
#' Ghorbani, M., Vafaei, N. and Myllymäki, M. (2025). A kernel-based test for the first-order separability
#' of spatio-temporal point processes, \emph{TEST}.
#'
#' @examples
#'
#' set.seed(123)
#' X <- cbind(runif(100), runif(100), runif(100, 0, 10))
#' perms <- block.permut(nblocks = 5, X = X, nperm = 10)
#' head(perms[[1]], 5)
#'
#' @export

block.permut <- function(nblocks, X, nperm = 1999) {
  # --- checks ---
  if (!is.numeric(nblocks) || length(nblocks) != 1 || nblocks < 2)
    stop("`nblocks` must be a single integer >= 2.")
  nblocks <- as.integer(nblocks)

  X <- as.matrix(X)
  if (ncol(X) < 3) stop("`X` must have at least 3 columns: x, y, time.")
  if (!is.numeric(X[,3])) stop("Time column X[,3] must be numeric.")
  if (anyNA(X[,3])) stop("Time column contains NA; please handle before calling.")
  if (!is.numeric(nperm) || length(nperm) != 1 || nperm < 1)
    stop("`nperm` must be a single integer >= 1.")
  nperm <- as.integer(nperm)

  # order by time
  ord <- order(X[,3])
  X <- X[ord, , drop = FALSE]

  N <- nrow(X)
  if (N < nblocks) stop("`nrow(X)` must be >= `nblocks`.")

  n <- N %/% nblocks
  size <- n * nblocks

  X_red <- X[seq_len(size), , drop = FALSE]
  tail_rows <- if (size < N) X[(size+1):N, , drop = FALSE] else NULL

  # build block indices (size x 1)
  # blocks are consecutive in time after ordering
  block_id <- rep(seq_len(nblocks), each = n)

  # store time values by block
  time_by_block <- split(X_red[,3], block_id)

  # number of all permutations (excluding identity)
  nall_minus1 <- factorial(nblocks) - 1

  # helper to generate ONE random non-identity permutation
  rand_perm <- function() {
    p <- sample.int(nblocks)
    while (all(p == seq_len(nblocks))) p <- sample.int(nblocks)
    p
  }

  # decide which permutations to use
  perms <- vector("list", length = min(nperm, nall_minus1))

  if (nall_minus1 <= nperm) {
    # enumerate all permutations if feasible (still can be huge; guard!)
    if (nblocks > 9) {
      stop("Enumerating all permutations is too large; reduce `nblocks` or set smaller `nperm`.")
    }
    # enumerate via recursion (small nblocks only)
    allp <- combinat::permn(seq_len(nblocks))  # if you keep combinat for small nblocks
    allp <- allp[!vapply(allp, function(p) all(p == seq_len(nblocks)), logical(1))]
    perms <- allp
  } else {
    # sample random unique permutations (approx; ensures no identity)
    # if you need strict uniqueness, add a hash set (costly). Typically unnecessary.
    for (i in seq_along(perms)) perms[[i]] <- rand_perm()
  }

  # generate permuted datasets
  out <- vector("list", length = length(perms))
  for (i in seq_along(perms)) {
    p <- perms[[i]]
    permuted_times <- unlist(time_by_block[p], use.names = FALSE)

    Xm <- X_red
    Xm[,3] <- permuted_times
    if (!is.null(tail_rows)) Xm <- rbind(Xm, tail_rows)

    out[[i]] <- Xm
  }

  out
}

