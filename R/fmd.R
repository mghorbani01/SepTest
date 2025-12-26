#' Spatio-temporal analysis of the 2001 FMD outbreak (illustrative example)
#'
#' This function test first-order separability for
#' the foot-and-mouth disease (FMD) outbreak data observed at the North Cumbria
#' study region provided by the \pkg{stpp} package.
#'
#' The example includes kernel-based estimates of the spatio-temporal intensity and related separability
#' diagnostics on a regular spatio-temporal grid for fmd data as well as the dHSIC test.
#'
#' @details
#' The datasets \code{fmd} and \code{northcumbria} are loaded from the
#' \pkg{stpp} package and are used purely for illustration.
#' They do not contain confidential information and have no scientific value.
#'
#' @return
#' An object of class \code{"test"} returned by
#' \code{\link{dHS.test}}.
#'
#' @references
#' Diggle, P., Rowlingson, B. and Su, T. (2005).
#' Point process methodology for on-line spatio-temporal disease surveillance.
#' \emph{Environmetrics}, 16, 423--434.
#'
#' @seealso
#' \code{\link[stpp]{fmd}},
#' \code{\link[stpp]{northcumbria}},
#' \code{\link{S.based.functions}},
#' \code{\link{dHS.test}}
#'
#' @examples
#'\dontrun{
#' if (requireNamespace("stpp", quietly = TRUE) ) {
#'
#'   library(stpp)
#'
#'   data("fmd", "northcumbria", package = "stpp")
#'
#'   s.region <- northcumbria / 1000
#'   t.region <- c(0, 200)
#'
#'   X <- as.3dpoints(
#'     fmd[, 1] / 1000,
#'     fmd[, 2] / 1000,
#'     fmd[, 3]
#'   )
#'
#'   ObsW <- spatstat.geom::owin(
#'     poly = list(
#'       x = s.region[, 1],
#'       y = s.region[, 2]
#'     )
#'   )
#'
#'   dHS.test(
#'     X,
#'     sim.procedure = "block_per",
#'     nblocks = 5L,
#'     nperm = 199L,
#'     nsim = 19L,
#'     bandwidth = NULL
#'   )
#'  }
#'}
#'
#' @export
run_fmd_example <- function() {
  ## implementation
}
