#' Plot a spatio-temporal point pattern
#'
#' Provides flexible visualization tools for spatio-temporal point patterns.
#' Depending on the chosen display type, the function produces one of three plots:
#' \itemize{
#'   \item A 3D scatterplot showing event locations in space–time (\code{"3D"});
#'   \item A 2D spatial projection of points in the spatial plane (\code{"space"});
#'   \item A temporal histogram with an overlaid kernel density curve (\code{"time"}).
#' }
#' The input dataset must contain three columns corresponding to the spatial (\code{x}, \code{y})
#' and temporal (\code{t}) coordinates of events.
#'
#' @param data A numeric matrix or data frame with at least three columns representing event coordinates.
#'   If \code{data} is a matrix, the first three columns are interpreted as \code{x}, \code{y}, and \code{t}.
#'   If \code{data} is a data frame, the function uses columns named \code{x}, \code{y}, \code{t} when present;
#'   otherwise it uses the first three columns.
#' @param type Character string specifying the type of visualization to produce:
#'   \code{"3D"} for a three-dimensional scatterplot;
#'   \code{"space"} for a 2D spatial plot of \code{x} versus \code{y};
#'   or \code{"time"} for a histogram of event times with a smoothed density overlay.
#' @param time_bins Integer specifying the number of bins in the histogram when \code{type = "time"}.
#'   Default is 30.
#' @param title Optional character string giving a plot title. If \code{NULL}, a default title is used.
#'
#' @details
#' The function serves as an exploratory tool for investigating the spatial,
#' temporal, or joint space–time structure of point pattern data.
#' Such visualization is often a first step before conducting statistical analyses
#' of separability or intensity modeling (see Ghorbani et al., 2021).
#'
#' \describe{
#'   \item{\strong{3D mode}}{Displays events in a 3D coordinate system using the
#'   \pkg{scatterplot3d} package, allowing a quick assessment of clustering
#'   and temporal trends in space--time.}
#'   \item{\strong{Spatial mode}}{Projects points onto the spatial plane (\code{x}--\code{y}),
#'   showing spatial structure independent of time.}
#'   \item{\strong{Temporal mode}}{Displays the marginal temporal distribution of
#'   events as a histogram with a kernel density overlay, facilitating the
#'   visual detection of temporal nonstationarity.}
#' }
#'
#' Visualization is particularly useful for validating the realism of simulated
#' data (e.g. from \code{\link{rstpoispp}} or \code{\link{Gauss.st.F}})
#' and for preliminary inspection prior to applying tests such as
#' \code{\link{chi2.test}}, \code{\link{global.envelope.test}}, or
#' \code{\link{dHS.test}}.
#'
#' @return
#' Produces a plot as a side effect. Nothing is returned.
#'
#' @note
#' The 3D visualization requires the \pkg{scatterplot3d} package.
#'
#' @importFrom graphics plot hist lines par
#' @importFrom stats density
#' @seealso
#' \code{\link{rstpoispp}},
#' \code{\link{st.intensity}},
#' \code{\link{plot_stlgcp}},
#' \code{\link{chi2.test}},
#' \code{\link{dHS.test}}
#'
#' @author
#' Nafiseh Vafaei \email{nafiseh.vafaei@slu.se}
#'
#' @references
#' Ghorbani M., Vafaei N., Dvořák J., Myllymäki M. (2021).
#' Testing the first-order separability hypothesis for spatio-temporal point patterns.
#' \emph{Computational Statistics & Data Analysis}, \bold{161}, 107245.
#'
#' @examples
#'
#' \donttest{
#' set.seed(123)
#' X <- cbind(runif(100), runif(100), runif(100, 0, 10))
#'
#' # Visualize point pattern in 3D space–time
#' plot_stpp(X, type = "3D")
#'
#' # View spatial projection
#' plot_stpp(X, type = "space")
#'
#' # Inspect temporal distribution
#' plot_stpp(X, type = "time", time_bins = 20)
#' }
#'
#' @export

 plot_stpp <- function(data,
                       type = c("3D", "space", "time"),
                       time_bins = 30,
                       title = NULL) {

   type <- match.arg(type)

   # ---- coerce + extract x,y,t safely ----
   if (is.matrix(data)) {
     if (ncol(data) < 3L) stop("`data` must have at least 3 columns (x, y, t).")
     data <- as.data.frame(data)
     x <- data[[1L]]; y <- data[[2L]]; tt <- data[[3L]]
   } else if (is.data.frame(data)) {
     if (all(c("x", "y", "t") %in% names(data))) {
       x <- data[["x"]]; y <- data[["y"]]; tt <- data[["t"]]
     } else {
       if (ncol(data) < 3L) stop("`data` must have at least 3 columns (x, y, t).")
       x <- data[[1L]]; y <- data[[2L]]; tt <- data[[3L]]
     }
   } else {
     stop("`data` must be a matrix or data.frame.")
   }

   if (!is.numeric(x) || !is.numeric(y) || !is.numeric(tt)) {
     stop("The first three columns (or columns x,y,t) must be numeric.")
   }
   if (length(x) == 0L) {
     warning("No points to plot.")
     return(invisible(NULL))
   }

   if (!is.numeric(time_bins) || length(time_bins) != 1L || !is.finite(time_bins) || time_bins < 1) {
     stop("`time_bins` must be a single integer >= 1.")
   }
   time_bins <- as.integer(time_bins)

   # default titles
   if (is.null(title)) {
     title <- switch(type,
                     "3D" = "Spatio-temporal point pattern",
                     "space" = "Spatial point pattern",
                     "time" = "Temporal distribution",
                     "Spatio-temporal point pattern"
     )
   }

   # ---- plotting ----
   if (type == "3D") {
     if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
       stop("Package 'scatterplot3d' is required for type = '3D'. Please install it.", call. = FALSE)
     }
     op <- graphics::par(mar = c(5, 5, 4, 2) + 0.1)
     on.exit(graphics::par(op), add = TRUE)

     scatterplot3d::scatterplot3d(
       x, y, tt,
       xlab = "x", ylab = "y", zlab = "t",
       pch = 20,
       angle = 40,
       box = FALSE,
       grid = TRUE,
       main = title
     )

   } else if (type == "space") {

     op <- graphics::par(mar = c(5, 5, 4, 2) + 0.1)
     on.exit(graphics::par(op), add = TRUE)

     graphics::plot(
       x, y,
       xlab = "x", ylab = "y",
       main = title,
       pch = 20
     )

   } else { # type == "time"

     op <- graphics::par(mar = c(5, 5, 4, 2) + 0.1)
     on.exit(graphics::par(op), add = TRUE)

     h <- graphics::hist(
       tt,
       breaks = time_bins,
       xlab = "Time",
       ylab = "Count",
       main = title,
       col = "gray",
       border = "white",
       freq = TRUE
     )

     # density overlay scaled to histogram counts
     d <- stats::density(tt)
     if (max(d$y) > 0) {
       scale_factor <- max(h$counts) / max(d$y)
       graphics::lines(d$x, d$y * scale_factor, col = "red", lwd = 2)
     }
   }

   invisible(NULL)
 }
