#' Plot spatio-temporal log-Gaussian Cox process (LGCP) realizations as time-sliced maps
#'
#' Produces a sequence of spatial raster maps displaying the evolution of a
#' simulated or fitted spatio-temporal log-Gaussian Cox process (LGCP) over time.
#' Each map shows the latent Gaussian random field (intensity surface) at a given
#' time slice, with the corresponding observed point events overlaid.
#' This visualization helps interpret temporal evolution and spatial clustering
#' patterns in simulated or fitted LGCP models.
#'
#' @param data A list containing components from a spatio-temporal LGCP model or simulation output:
#' \describe{
#'   \item{\code{RF}}{A list describing the latent Gaussian random field, with elements
#'         \code{xcoord}, \code{ycoord}, \code{tcoord}, and \code{Z}, where \code{Z}
#'         is a 3D array of intensity (or log-intensity) surfaces over space and time.}
#'   \item{\code{st.lgcp}}{A data frame with columns \code{x}, \code{y}, and \code{t}
#'         giving the spatial and temporal coordinates of observed events.}
#' }
#'
#' @details
#' The function plots up to 10 evenly spaced time slices from the latent intensity field
#' and overlays the corresponding point events accumulated up to each time point.
#' Each panel represents a spatial realization at a fixed time \eqn{t_k}, providing
#' a visual summary of the dynamic spatio-temporal structure of the LGCP.
#'
#' This approach follows the visualization principles used in Ghorbani et al. (2021, 2025),
#' where LGCPs are employed to assess the behavior of separability diagnostics and
#' kernel-based tests in complex, non-separable point process settings.
#' The raster intensity surfaces illustrate latent heterogeneity, while the overlaid
#' points display observed event clustering relative to the underlying field.
#'
#' Time slices are selected at evenly spaced quantiles of the temporal domain, and
#' each plot includes:
#' \itemize{
#'   \item A raster map of the latent intensity surface at time \eqn{t_k};
#'   \item White contour lines showing equal-intensity regions;
#'   \item Overlaid black points representing events observed up to \eqn{t_k}.
#' }
#'
#' The resulting plots are arranged into a single grid layout using the
#' \pkg{patchwork} package for ease of comparison.
#'
#' @return
#' A combined \code{ggplot} object displaying up to ten spatial raster maps arranged
#' in a grid layout (by default, two rows and up to five columns).
#' Each panel corresponds to one time slice. The function returns the combined plot object.
#'
#' @importFrom ggplot2 aes geom_raster geom_contour geom_point scale_fill_viridis_c
#' @importFrom ggplot2 labs coord_fixed theme_minimal theme
#' @importFrom reshape2 melt
#' @importFrom patchwork wrap_plots
#' @importFrom graphics par
#' @import ggplot2
#'
#' @seealso
#' \code{\link{Gauss.st.F}} for simulating spatio-temporal Gaussian random fields;
#' \code{\link{get.lambda.function}} and \code{\link{rstpoispp}}
#' for generating intensity-based spatio-temporal point processes.
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
#' \donttest{
#' # Example: visualize a spatio-temporal LGCP simulation
#' out <- rstLGCPP(xlim = c(0,1),
#'                 ylim = c(0,1),
#'                 tlim = c(0,1),
#'                 grid = c(15,15,10))
#' plot_stlgcp(data = out)
#'}
#' @export
plot_stlgcp <- function(data) {

  # ---- dependency checks (CRAN-friendly) ----
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_stlgcp(). Please install it.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for plot_stlgcp(). Please install it.", call. = FALSE)
  }

  # ---- validate structure ----
  if (!is.list(data) || is.null(data$RF) || is.null(data$st.lgcp)) {
    stop("`data` must be a list with components `RF` and `st.lgcp`.", call. = FALSE)
  }

  g <- data$RF
  pts <- data$st.lgcp

  if (!is.list(g) || any(vapply(c("xcoord", "ycoord", "tcoord", "Z"), function(nm) is.null(g[[nm]]), logical(1)))) {
    stop("`data$RF` must contain `xcoord`, `ycoord`, `tcoord`, and `Z`.", call. = FALSE)
  }
  if (!is.data.frame(pts) || !all(c("x", "y", "t") %in% names(pts))) {
    stop("`data$st.lgcp` must be a data.frame with columns `x`, `y`, `t`.", call. = FALSE)
  }

  x <- g$xcoord
  y <- g$ycoord
  tt <- g$tcoord
  Z <- g$Z

  if (!is.numeric(x) || !is.numeric(y) || !is.numeric(tt)) stop("`xcoord`, `ycoord`, and `tcoord` must be numeric.")
  if (!is.array(Z) || length(dim(Z)) != 3L) stop("`RF$Z` must be a 3D array.")
  if (dim(Z)[1L] != length(x) || dim(Z)[2L] != length(y) || dim(Z)[3L] != length(tt)) {
    stop("Dimensions of `RF$Z` must match lengths of `xcoord`, `ycoord`, and `tcoord`.")
  }

  nt <- length(tt)
  if (nt < 1L) stop("`RF$tcoord` must have positive length.")

  # ---- choose up to 10 time indices, evenly spaced, unique ----
  k_to_plot <- unique(round(seq.int(1, nt, length.out = min(10L, nt))))
  k_to_plot <- k_to_plot[k_to_plot >= 1L & k_to_plot <= nt]

  plots <- list()

  for (k in k_to_plot) {
    sel <- pts$t <= tt[k]
    if (!any(sel)) next

    Z_slice <- Z[, , k, drop = TRUE]

    # base-R melt: as.data.frame(as.table()) avoids reshape2 dependency
    df <- as.data.frame(as.table(Z_slice), responseName = "z")
    # table indices are factors; convert to integer
    i <- as.integer(df$Var1)
    j <- as.integer(df$Var2)

    df$x <- x[i]
    df$y <- y[j]
    df <- df[is.finite(df$z), , drop = FALSE]

    pts_k <- data.frame(x = pts$x[sel], y = pts$y[sel])

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = z), interpolate = TRUE) +
      ggplot2::geom_contour(ggplot2::aes(z = z), color = "white", bins = 10, linewidth = 0.3) +
      ggplot2::geom_point(data = pts_k, ggplot2::aes(x = x, y = y),
                          inherit.aes = FALSE, color = "black", size = 0.8) +
      ggplot2::scale_fill_viridis_c(option = "magma") +
      ggplot2::labs(title = sprintf("t = %.2f", tt[k]), x = "x", y = "y") +
      ggplot2::coord_fixed() +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "none")

    plots[[length(plots) + 1L]] <- p
  }

  if (length(plots) == 0L) {
    stop("No time slices contained any points (check `st.lgcp$t` vs `RF$tcoord`).", call. = FALSE)
  }

  patchwork::wrap_plots(plots, ncol = min(5L, length(plots)))
}
