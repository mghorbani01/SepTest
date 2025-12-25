#' Simulate an inhomogeneous spatio-temporal Poisson point process
#'
#' Generates a realization of an inhomogeneous Poisson point process (STPP)
#' in space and time using the standard thinning method.
#' The user provides an intensity function \eqn{\lambda(u, t)} and an upper
#' bound \eqn{L_{\max}} on its value over the observation window.
#' The algorithm first samples candidate events uniformly over space and time
#' and then retains each candidate with probability proportional to its
#' normalized intensity \eqn{\lambda(u,t)/L_{\max}}.
#'
#' @param lambda A function of the form \code{lambda(u, t)} that returns the intensity value at coordinates \code{(x, y, t)}.
#' @param Lmax A numeric value giving the known or estimated maximum of the intensity function \code{lambda} over the spatial and temporal window. Used for thinning.
#' @param s.region A matrix with two columns giving the polygonal spatial window. Each row is a vertex of the polygon. Default is the unit square.
#' @param t.region A numeric vector of length 2 giving the temporal observation window. Default is \code{c(0,1)}.
#'
#'
#' @details
#' The method implements the classical \emph{thinning algorithm} for simulating
#' inhomogeneous Poisson processes:
#'
#' \enumerate{
#'   \item Draw \eqn{N^* \sim \mathrm{Poisson}(L_{\max} \, |W| \, |T|)}, where
#'   \eqn{|W|} and \eqn{|T|} denote the spatial and temporal window measures.
#'   \item Generate \eqn{N^*} candidate points uniformly over
#'   \eqn{W \times T}.
#'   \item Retain each point \eqn{(u_i, t_i)} independently with probability
#'   \eqn{p_i = \lambda(u_i, t_i) / L_{\max}}.
#' }
#'
#' The result is a realization of an inhomogeneous STPP with intensity function
#' \eqn{\lambda(u, t)}.
#'
#'
#' This simulator underpins the spatio-temporal framework introduced in
#' Ghorbani et al. (2021, 2025) for studying \emph{first-order separability}.
#' By selecting appropriate intensity functions (see \code{\link{get.lambda.function}}),
#' users can generate fully separable, partially separable, or non-separable
#' spatio-temporal patterns, enabling direct evaluation of separability tests such as
#' \code{\link{chi2.test}}, \code{\link{global.envelope.test}}, or
#' \code{\link{dHS.test}}.
#'
#' @return A numeric matrix with three columns (\code{x}, \code{y}, \code{t}) representing the retained points from the inhomogeneous Poisson process.
#'
#' @note
#' The intensity function \eqn{\lambda(u, t)} should return non-negative
#' numeric values and be bounded above by \code{Lmax} across the observation domain.
#'
#' @importFrom stats rpois runif
#' @importFrom splancs as.points
#' @seealso
#' \code{\link{get.lambda.function}} to construct spatio-temporal intensity models;
#' \code{\link{get.lambda.max}} to compute intensity maxima;
#' \code{\link{estimate.st.intensity}} for intensity estimation;
#' \code{\link{plot_stpp}} for visualization.
#' @export
#'
#'
#' @author
#' Mohammad Ghorbani \email{mohammad.ghorbani@slu.se}\cr
#' Nafiseh Vafaei \email{nafiseh.vafaei@ltu.se}
#'
#' @references
#' Ghorbani M., Vafaei N., Dvořák J., Myllymäki M. (2021).
#' Testing the first-order separability hypothesis for spatio-temporal point patterns.
#' \emph{Computational Statistics & Data Analysis}, \bold{161}, 107245.
#'
#' Ghorbani, M., Vafaei, N. and Myllymäki, M. (2025). A kernel-based test for the first-order separability
#' of spatio-temporal point processes, \emph{TEST} .
#'
#' @examples
#'
#' \donttest{
#' # Example 1: Simulate a separable spatio-temporal Poisson process
#' lambda <- get.lambda.function(N = 200, g = 50, model = 1)
#' Lmax   <- get.lambda.max(N = 200, g = 50, model = 1)
#' X <- rstpoispp(lambda, Lmax)
#' head(X)
#'
#' # Example 2: Non-separable model (Model 4)
#' lambda <- get.lambda.function(N = 200, g = 50, model = 4)
#' Lmax   <- get.lambda.max(N = 200, g = 50, model = 4)
#' sim_data <- rstpoispp(lambda, Lmax)
#'
#' # Spatial projection of simulated events
#' plot(sim_data[, 1:2], asp = 1, main = "Spatial Projection of Simulated stPP")
#' # Example 3: 3D visualization using plot_ST_pp()
#' plot_stpp(X, type = "3D", title="Realisation of a stPP")
#' }

rstpoispp <- function(lambda, Lmax,
                        s.region=splancs::as.points(c(0,1,1,0),c(0,0,1,1)),
                        t.region=c(0,1)) {
  n <- rpois(1, Lmax)
  x <- runif(n, min(s.region[,1]), max(s.region[,1]) )
  y <- runif(n, min(s.region[,2]), max(s.region[,2]) )
  t <- runif(n, t.region[1], t.region[2])
  prob <- lambda(x,y,t)/Lmax
  u <- runif(n)
  retain <- (u <= prob)
  xyt <- cbind(x,y,t)
  out <- xyt[retain,]
  return(out)
}
