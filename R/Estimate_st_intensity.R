#' Kernel Estimation of the Spatio-Temporal Intensity Function and Its Components,
#' and Test Statistics for First-Order Separability
#'
#' It returns estimates of the full spatio-temporal intensity together with its
#' spatial and temporal marginal components.
#' The output also includes the test statistic
#' \eqn{S(u,t)} and its spatial and temporal marginal profiles
#' \eqn{S_{\mathrm{space}}(u)} and \eqn{S_{\mathrm{time}}(t)}
#' (equations (11)–(13) in Ghorbani et al., 2021),
#' as well as several deviation tests (e.g., equation (15) in Ghorbani et al., 2021)
#' that quantify departures from first-order separability.
#'
#' @param X A numeric matrix with three columns giving the event coordinates
#' (\code{x}, \code{y}, and \code{t}).
#' @param s.region A numeric matrix with two columns defining the polygonal
#' boundary of the spatial observation region.
#' @param t.region A numeric vector of length 2 specifying the temporal observation window.
#' @param at Character string; either \code{"pixels"} to estimate intensity on a
#' spatio-temporal grid or \code{"points"} to compute the estimates at observed event locations.
#' @param n.grid A numeric vector of length 3 giving the number of grid cells in
#' the \code{x}, \code{y}, and \code{t} dimensions. Required when \code{at = "pixels"}.
#'
#' @details
#' The estimation follows the kernel framework of
#' Ghorbani et al. (2021).  A Gaussian kernel is applied in each
#' spatial and temporal dimension.  Spatial bandwidths are estimated using
#' Diggle’s method, and the temporal bandwidth is obtained by the
#' Sheather–Jones direct plug-in (SJ-DPI) approach
#' (see \code{\link{calc.bandwidths.and.edgecorr}} for implementation details).
#' Edge corrections are computed analytically using Gaussian tail probabilities
#'  (see \code{\link{calc.bandwidths.and.edgecorr}} for implementation details).
#'
#'
#' The non-separable intensity estimate is
#' \deqn{
#' \hat{\rho}(u,t) =
#' \sum_{i=1}^n
#'   \frac{K_\epsilon^2\!\left(u-u_i\right)}{C_{W,\epsilon}(u_i)}
#'   \frac{K_\delta^1\!\left(t-t_i\right)}{C_{T,\delta}(t_i)},
#' }
#' where \eqn{k_b(v) = k(v/b)/b^d},  and \eqn{K^1} and \eqn{K^2} are Gaussian kernels with spatial and temporal
#' bandwidths \eqn{\epsilon} and \eqn{\delta}. \eqn{C_{W,\epsilon}(u_i)} and \eqn{C_{T,\delta}(t_i)} are spatial and temporal edge corrections, respectively.
#' For estimate of the separable counterparts and more details see equations (3)-(5) in Ghorbani et al., 2021).
#'
#' The separability test is:
#' \deqn{S(u, t) = \frac{\hat \rho(u, t)}{\hat\rho_{\mathrm{space}}(u)\hat\rho_{\mathrm{time}}(t)/n},}
#' where \eqn{n} is the number of observed points.
#'
#' Several deviation statistics are provided to quantify departures from separability
#' \itemize{
#'   \item \code{deviation.t1}: Integral of absolute deviations \eqn{\int_{W\times T}|\hat\rho(u,t) - \hat\rho_{\mathrm{sep}}(u,t)|\mathrm{d}\,u\mathrm{d}\,t}.
#'   \item \code{deviation.t2}: Integral of absolute deviation of inverse intensities  \eqn{\int_{W\times T}|1/\hat\rho(u,t) - 1/\hat\rho_{\mathrm{sep}}(u,t)|\mathrm{d}\,u\mathrm{d}\,t}.
#'   \item \code{deviation.t3}: Sum of log-ratio deviations \eqn{\sum_i\big(\log(\hat\rho(u_i,t_i)- \log(\hat\rho_{\mathrm{sep}}(u_i,t_i))\big)}.
#'   \item \code{deviation.t4}: Integral of the S-function.
#' }
#'
#' When \code{at = "points"}, intensities and diagnostics are evaluated at
#' the observed event locations.
#'
#' @return A list containing:
#' \describe{
#'   \item{x, y, t}{Grid vectors for x, y, and t (returned only when \code{at = "pixels"}).}
#'   \item{epsilon}{Estimated spatial bandwidth used in Gaussian kernels.}
#'   \item{delta}{Estimated temporal bandwidth.}
#'   \item{SPat.intens}{Estimated spatial intensity surface (pixels only).}
#'   \item{TeM.intens}{Estimated temporal intensity profile (pixels only).}
#'   \item{sep.intens}{Separable spatio-temporal intensity estimate.}
#'   \item{nonsep.intens}{Non-separable spatio-temporal intensity estimate.}
#'   \item{S.fun}{Test function \eqn{S(u,t)} as the ratio of non-separable to separable intensity estimates.}
#'   \item{S.space}{Marginal sum of \eqn{S(u,t)} over time (space profile).}
#'   \item{S.time}{Marginal sum of \eqn{S(u,t)} over space (time profile).}
#'   \item{deviation.t1}{Deviation statistic: Integral of absolute deviations.}
#'   \item{deviation.t2}{Deviation statistic: absolute deviation of inverse intensities.}
#'   \item{deviation.t3}{Deviation statistic: sum of log-ratio deviations.}
#'   \item{deviation.t4}{Total integral of the S-function.}
#' }
#'
#' @seealso
#' \code{\link{calc.bandwidths.and.edgecorr}},
#' \code{\link{global.envelope.test}},
#' \code{\link{chi2.test}},
#' \code{\link{dHS.test}}
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
#' @examples
#' \donttest{
#' set.seed(123)
#' X <- cbind(runif(100), runif(100), runif(100, 0, 10))
#' s.region <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), ncol = 2, byrow = TRUE)
#' t.region <- c(0, 10)
#' result <- st.intensity(X, s.region, t.region, at = "pixels", n.grid = c(64, 64, 32))
#' str(result$S.fun)
#' image(result$S.fun[,,5], main = "S-function slice at t[5]",
#'        col = topo.colors(50))
#' contour(result$S.fun[,,5], add=TRUE)
#' str(result[c("deviation.t1","deviation.t3","deviation.t4")])
#' }
#'
#' @export
st.intensity <-function(X, s.region, t.region, at=c("pixels","points"),
                        n.grid = c(25,25,20)) {
  at <- match.arg(at)  # This ensures 'at' is either "pixels" or "points" and length 1
  np <- length(X[,1]) # Number of points
  edge <- calc.bandwidths.and.edgecorr(X, s.region, t.region, n.grid)

  switch(at,
         pixels = {
           g.s <-(max(s.region[,1])-min(s.region[,1]))/n.grid[1]*(max(s.region[,2])-min(s.region[,2]))/n.grid[2] # grid size in space
           t.s <-((t.region[2]-t.region[1])/n.grid[3]) # partiton size in time axis.
           # Sizes of the grid cells
           dx <- (max(s.region[,1])-min(s.region[,1]))/n.grid[1]
           dy <- (max(s.region[,2])-min(s.region[,2]))/n.grid[2]
           dt <- (max(t.region)-min(t.region))/n.grid[3]
           # Centers of the pixels
           gx <- seq.int(min(s.region[,1]) + dx/2, max(s.region[,1]) - dx/2, length.out = n.grid[1]) # the center of the grids at x axis
           gy <- seq.int(min(s.region[,2]) + dy/2, max(s.region[,2]) - dy/2, length.out = n.grid[2]) # the center of the grids at y axis
           gt <- seq.int(min(t.region) + dt/2, max(t.region) - dt/2, length.out = n.grid[3])         # the center of the grids at t axis

           ax <- outer(gx, X[,1], dnorm, sd = edge$bw[1L])/matrix(rep(edge$space, n.grid[1]), ncol=np, nrow=n.grid[1], byrow=T) # the contribution of each point
           # in intensity estimation
           ay <- outer(gy, X[,2], dnorm, sd = edge$bw[2L])
           at <- outer(gt, X[,3], dnorm, sd = edge$bw[3L])/matrix(rep(edge$time, n.grid[3]), ncol=np, nrow=n.grid[3], byrow=T)

           # spatial intensity estimation
           Spat.intens <- tcrossprod(ax,ay) #t(ax)%*%ay
           Spat.intens <- pmax(Spat.intens, 2^-1024)
           # temporal intensity estimation
           temp.intens <- apply(at, 1, sum)
           # seprable spatio-temporal intensity estimation
           SEP.intensity <- (Spat.intens%o%temp.intens)/np
           SEP.intensity <- pmax(SEP.intensity, 2^-1024)
           # non-seprable spatio-temporal intensity estimation
           non.sep <- array(0, n.grid)
           t.ay <- t(ay)
           for (k in 1:(n.grid[3])) {
             t.ay.tk <- t.ay * at[k, ]
             non.sep[, ,k ] <- ax %*% t.ay.tk
           }
           ## test statistics
           S.fun <- non.sep/SEP.intensity
           S.space <- apply(S.fun, 1:2, sum)*t.s     # space
           S.time <- apply(S.fun, 3, sum)*g.s
           # Deviation tests
           Sd1 <- sum(abs(non.sep-SEP.intensity))*t.s*g.s
           Sd2 <- sum(abs(1/non.sep - 1/SEP.intensity))*t.s*g.s
           Sd3 <- sum(log(non.sep)-log(SEP.intensity))*t.s*g.s
           Sd4 <- sum(S.fun)*t.s*g.s
         },
         points={
           ax <- outer(X[,1], X[,1], dnorm, sd = edge$bw[1L])/matrix(rep(edge$space, np), ncol=np,  byrow=T)
           ay <- outer(X[,2], X[,2], dnorm, sd = edge$bw[2L])
           at <- outer(X[,3], X[,3], dnorm, sd = edge$bw[3L])/matrix(rep(edge$time, np), ncol=np, byrow=T)

           #Spat.intens  <- spatstat::density.ppp(SPP, at="points", sigma=spatstat::bw.diggle(SPP), edge=TRUE, diggle=TRUE,leaveoneout = T)
           Spat.intens <- ax*t(ay)
           diag(Spat.intens) <- 0
           Spat.intens <- apply(Spat.intens, 1,sum)
           diag(at) <- 0
           temp.intens <- apply(at, 1, sum)

           ## separable intensity function
           SEP.intensity <- (Spat.intens*temp.intens)/np

           ### non=separable intensity function
           non.sep <- ax * t(ay) * at
           diag(non.sep)<-0
           non.sep <- apply(non.sep,1,sum)
           S.fun <- non.sep/SEP.intensity

           ####deviation tests at points
           # Deviation test
           Sd1 <-  sum(abs(non.sep-SEP.intensity))
           Sd2 <-  sum(abs(1/non.sep - 1/SEP.intensity))
           Sd3 <-  sum(log(non.sep)-log(SEP.intensity))
           Sd4 <- sum(S.fun)
           gx <-NULL
           gy <-NULL
           gt <-NULL
           S.space <- NULL
           S.time  <- NULL
         })

  return(list(x=gx, y=gy, t=gt, epsilon=edge$bw[1], delta=edge$bw[3],
              SPat.intens=Spat.intens, TeM.intens=temp.intens,
              nonsep.intens=non.sep, sep.intens=SEP.intensity,
              S.fun=S.fun,
              S.space=S.space, S.time=S.time,
              deviation.t1=Sd1,
              deviation.t2=Sd2,
              deviation.t3=Sd3,
              deviation.t4=Sd4))
}
