#' Plot spatio-temporal determinantal point process (DPP) realizations
#'
#' Produces diagnostic plots for spatio-temporal determinantal point process
#' (DPP) simulations or fitted models. The function formats the plot title
#' based on the spatial and temporal interaction parameters \eqn{\alpha_s} and
#' \eqn{\alpha_t}, automatically displaying either the parameters themselves or
#' their reciprocals when greater than 1.
#'
#' @param data A spatio-temporal point pattern object suitable for
#'   \code{plot_stpp()}, typically from the \pkg{stpp} or related packages.
#' @param type Character string specifying the type of plot to produce.
#'   This is passed directly to \code{plot_stpp()}. Typical values are \code{"3D"}, \code{"space"}, and \code{"time"}.
#' @param alpha_s Numeric scalar (> 0). Spatial interaction parameter of the DPP model.
#' @param alpha_t Numeric scalar (> 0). Temporal interaction parameter of the DPP model.
#'
#' @return No return value. The function is called for its side effect of
#'   producing a diagnostic plot.
#'
#' @details
#' If \eqn{\alpha_s > 1} and \eqn{\alpha_t > 1}, the title displays
#' \eqn{\alpha_s^{-1}} and \eqn{\alpha_t^{-1}}, which correspond to interaction
#' ranges. Otherwise, the parameters are shown directly.
#'
#' The function then calls \code{plot.ST.pp()} to produce the actual plot.
#'
#' @examples
#' \donttest{
#' # Simulate a stationary separable Matérn ST-DPP
#' sim <- rstDPP(
#'   mode     = "stationary",
#'   model    = "S",
#'   spectral = "matern",
#'   alpha_s  = 10,
#'   alpha_t  = 4.7,
#'   nu       = 2,
#'   eps      = 1,
#'   lambda_max = 70,
#'   grid_size  = 1.5
#' )
#'  plot_stDPP(sim, type = "3D",  alpha_s = 10, alpha_t = 4.7)
#' }
#'
#' @export
plot_stDPP <- function(data,
                        type = c("3D", "space", "time"),
                        alpha_s,
                        alpha_t) {

  type <- match.arg(type)

  if (!is.numeric(alpha_s) || length(alpha_s) != 1L || !is.finite(alpha_s) || alpha_s <= 0) {
    stop("`alpha_s` must be a single finite numeric value > 0.")
  }
  if (!is.numeric(alpha_t) || length(alpha_t) != 1L || !is.finite(alpha_t) || alpha_t <= 0) {
    stop("`alpha_t` must be a single finite numeric value > 0.")
  }

  if (alpha_s > 1 && alpha_t > 1) {
    main_title <- bquote(
      alpha[s]^{-1} == .(round(1 / alpha_s, 2)) ~ "," ~
        alpha[t]^{-1} == .(round(1 / alpha_t, 2))
    )
  } else {
    main_title <- bquote(
      alpha[s] == .(round(alpha_s, 2)) ~ "," ~
        alpha[t] == .(round(alpha_t, 2))
    )
  }

  plot_stpp(
    data  = data,
    type  = type,
    title = main_title
  )

  invisible(NULL)
}

