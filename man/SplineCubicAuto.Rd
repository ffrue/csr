% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CubicSplineRegression.R
\name{SplineCubicAuto}
\alias{SplineCubicAuto}
\title{Fit Upper Frontier Cubic Spline with Optimal Number of Knots}
\usage{
SplineCubicAuto(
  xtab,
  ytab,
  xgrid,
  cv,
  ic = "BIC",
  kmax = 5,
  iters = 10000,
  eps = 1e-06
)
}
\arguments{
\item{xtab}{numeric vector of x-indices of data points}

\item{ytab}{numeric vector of y-indices of data points}

\item{xgrid}{numeric vector of x-values where function should be evaluated}

\item{cv}{integer. indicates constraints for fitting: -1 (unconstrained), 0 (monotonically increasing), 1 (monotonically increasing + concave)}

\item{ic}{string. information criterion: either "AIC" or "BIC"}

\item{kmax}{integer. highest number of knots considered by information criterion}

\item{iters}{integer. iterations used for optimization for the CVXR solver}

\item{eps}{integer. stopping criterion (absolute some of deviations from constraints)}
}
\value{
A list with:
\describe{
\item{\code{fit}}{Numeric vector of fitted y-values at the grid points.}
\item{\code{k}}{The selected number of knots.}
}
}
\description{
fits an (un)constrained cubic spline as upper frontier, whose number of knots is chosen optimally by an information criterion.
}
\details{
This function fits cubic splines with different numbers of knots (from 1 to the kmax chosen)
and then selects the one with the lowest AIC or BIC value. The information criteria are calculated as
\deqn{AIC = \log(\sum|y-i-\hat{phi}(x_i)|) + 2 \cdot (k+4) / n} and
\deqn{BIC = \log(\sum|y-i-\hat{phi}(x_i)|) + \log(n) \cdot (k+4) / n}
}
\examples{
xtab = atp_2017$ranking_points
ytab = atp_2017$wins
xgrid = seq(min(xtab),max(xtab),length.out=101)
fit_m <- SplineCubicAuto(xtab=xtab, ytab=ytab, xgrid=xgrid,
                cv=0, ic="BIC", kmax=3, iters=1000000, eps=1e-6)
fit_mc <- SplineCubicAuto(xtab=xtab, ytab=ytab, xgrid=xgrid,
                cv=1, ic="BIC", kmax=3, iters=1000000, eps=1e-6)
plot(xtab,ytab, xlab="ranking points", ylab="wins", main="ATP results in 2017 per player")
lines(xgrid,fit_m$fit,col="red")
lines(xgrid,fit_mc$fit,col="blue")
legend("bottomright", col = c("blue", "red"),
        legend = c("monontone+concave", "monotone"), lwd = 2)
}
