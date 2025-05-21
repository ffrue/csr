#' Find FDH Frontier Points
#'
#' finds the x-indices of points on the FDH-frontier (Free Disposal Hull)
#'
#' @param xtab numeric vector of x-indices of data points
#' @param ytab numeric vector of y-indices of data points
#' @return numeric vector with x-indices of points on the FDH-frontier
#' @details
#' The free disposal hull method was developed in Deprins, Simar and Tulkens (1984). It is the lowest stair-case curve on or above all data points and by design monotonically increasing.
#' Formally, it is defined as \deqn{\phi(x)_n = \max \{ y_i, i:x_i \leq x \} }
#' The function returns the indices of all data points on the FDH (where \eqn{y_i = \phi(x_i)}).
#' By linearly connecting these points one arrives at the linearized FDH curve.
#' @export
#'
#' @examples
#' xtab=atp_2017$ranking_points
#' ytab=atp_2017$wins
#' plot(xtab,ytab, xlab="ranking points", ylab="wins", main="ATP results in 2017 per player")
#' lines(xtab[fdh_index(xtab,ytab)],
#'         ytab[fdh_index(xtab,ytab)],
#'           type = "l", col = "red", lty = 2)
fdh_index <- function(xtab, ytab) {
  fidx <- c()

  for (i in 1:length(xtab)) {
    xcan <- xtab[i]
    feasible_indices <- which(xtab <= xcan)
    max_y <- max(ytab[feasible_indices])
    is_fdh <- ytab[i] == max_y

    if (is_fdh) {
      fidx <- c(fidx, i)
    }
  }

  fidx <- sort(fidx)
  return(fidx)
}

#' Find FDH Lower Frontier Points
#'
#' finds the x-indices of points on the lower FDH-frontier
#'
#' @param xtab numeric vector of x-indices of data points
#' @param ytab numeric vector of y-indices of data points
#' @return numeric vector with x-indices of points on the lower FDH-frontier
#' @details
#' The free disposal hull method was developed in Deprins, Simar and Tulkens (1984). It is the lowest stair-case curve on or below all data points and by design monotonically increasing.
#' Formally, it is defined as \deqn{\phi(x)_n = \min \{ y_i, i:x_i \geq x \} }
#' The function returns the indices of all data points on the FDH (where \eqn{y_i = \phi(x_i)}).
#' By linearly connecting these points one arrives at the linearized FDH curve.
#' @export
#'
#' @examples
#' xtab=atp_2017$ranking_points
#' ytab=atp_2017$wins
#' plot(xtab,ytab, xlab="ranking points", ylab="wins", main="ATP results in 2017 per player")
#' lines(xtab[fdh_lower_index(xtab,ytab)],
#'         ytab[fdh_lower_index(xtab,ytab)],
#'           type = "l", col = "red", lty = 2)
fdh_lower_index <- function(xtab, ytab) {
  fidx <- c()

  for (i in seq_along(xtab)) {
    xcan <- xtab[i]
    indices <- which(xtab >= xcan)
    is_fdh <- ytab[i] == min(ytab[indices])

    if (is_fdh) {
      fidx <- c(fidx, i)
    }
  }

  fidx <- sort(fidx)
  return(fidx)
}



#' Find DEA Frontier Points
#'
#' finds the x-indices of points on the DEA-frontier (with varying returns to scale)
#'
#' @param xtab numeric vector of x-indices of data points
#' @param ytab numeric vector of y-indices of data points
#' @return numeric vector with x-indices of points on the DEA-frontier
#' @details
#' The Data Envelopment Analysis (DEA) frontier under varying returns to scale (VRS) envelops all data points from above, forming a piecewise linear and monotonically increasing frontier.
#' This implementation follows the standard VRS approach as introduced by Banker, Charnes, and Cooper (1984).
#' The function returns the indices of all data points that lie on this frontier.
#' By linearly connecting these points, one obtains the VRS DEA frontier.
#' @export
#'
#' @examples
#' xtab=atp_2017$ranking_points
#' ytab=atp_2017$wins
#' plot(xtab,ytab, xlab="ranking points", ylab="wins", main="ATP results in 2017 per player")
#' lines(xtab[dea_index(xtab,ytab)],
#'         ytab[dea_index(xtab,ytab)],
#'           type = "l", col = "red", lty = 2)
dea_index <- function(xtab, ytab) {
  fidx <- fdh_index(xtab, ytab)

  xfdh <- xtab[fidx]
  yfdh <- ytab[fidx]

  ridx <- order(xfdh)
  xsfdh <- xfdh[ridx]
  ysfdh <- yfdh[ridx]

  didx <- c(1)
  xcur <- xsfdh
  ycur <- ysfdh

  while (length(xcur) >= 2) {
    tlength <- length(xcur)
    slopes <- (ycur[2:tlength] - ycur[1]) / (xcur[2:tlength] - xcur[1])
    b <- which.max(slopes)

    didx <- c(didx, didx[length(didx)] + b)

    xcur <- xcur[(b + 1):tlength]
    ycur <- ycur[(b + 1):tlength]
  }

  didx <- fidx[ridx[didx]]
  didx <- sort(didx)

  return(didx)
}

#' Find DEA Lower Frontier Points
#'
#' finds the x-indices of points on the lower DEA-frontier (with varying returns to scale)
#'
#' @param xtab numeric vector of x-indices of data points
#' @param ytab numeric vector of y-indices of data points
#' @return numeric vector with x-indices of points on the lower DEA-frontier
#' The Data Envelopment Analysis (DEA) frontier under varying returns to scale (VRS) envelops all data points from below, forming a piecewise linear and monotonically increasing frontier.
#' This implementation follows the standard VRS approach as introduced by Banker, Charnes, and Cooper (1984).
#' The function returns the indices of all data points that lie on this frontier.
#' By linearly connecting these points, one obtains the lower VRS DEA frontier.
#' @export
#'
#' @examples
#' xtab=atp_2017$ranking_points
#' ytab=atp_2017$wins
#' plot(xtab,ytab, xlab="ranking points", ylab="wins", main="ATP results in 2017 per player")
#' lines(xtab[dea_lower_index(xtab,ytab)],
#'         ytab[dea_lower_index(xtab,ytab)],
#'           type = "l", col = "red", lty = 2)
dea_lower_index <- function(xtab, ytab) {
  fidx <- fdh_lower_index(xtab, ytab)

  xfdh <- xtab[fidx]
  yfdh <- ytab[fidx]

  ord <- order(xfdh)
  xsfdh <- xfdh[ord]
  ysfdh <- yfdh[ord]

  didx <- 1
  xcur <- xsfdh
  ycur <- ysfdh

  while (length(xcur) >= 2) {
    tlength <- length(xcur)

    # Compute slopes
    slopes <- (ycur[2:tlength] - ycur[1]) / (xcur[2:tlength] - xcur[1])
    b <- which.min(slopes)

    didx <- c(didx, didx[length(didx)] + b)

    xcur <- xcur[(b + 1):tlength]
    ycur <- ycur[(b + 1):tlength]
  }

  final_didx <- fidx[ord[didx]]
  final_didx <- sort(final_didx)

  return(final_didx)
}


#' Fit Upper Frontier Cubic Spline with Predefined Number of Knots
#'
#' fits an (un)constrained cubic spline as upper frontier with predefined number of knots
#'
#' @param xtab numeric vector of x-indices of data points
#' @param ytab numeric vector of y-indices of data points
#' @param xgrid numeric vector of x-values where function should be evaluated
#' @param kn integer. number of knots for spline fitting
#' @param cv integer. indicates constraints for fitting: -1 (unconstrained), 0 (monotonically increasing), 1 (monotonically increasing + concave)
#' @param iters integer. iterations used for optimization for the CVXR solver
#' @param eps integer. stopping criterion (absolute some of deviations from constraints)
#' @return numeric vector with y-values of fitted cubic spline at grid points
#' @details The upper boundary curve fitted by this function is a cubic spline with selected number of knots \eqn{t_j} (location chosen at quantiles of $x$)
#' and the solution to the optimization problem \deqn{\min_\alpha \int \pi(x)^T \alpha dx = \sum_{j=0}^{N-1} \alpha_j \int \pi_j(x) dx}
#' subject to the envelopment constraint \eqn{\pi(x_i)^T \alpha \geq y_i \forall i=1,...,n}.
#' Intuitively, the boundary is chosen to minimize the area under the fitted curve while also being weakly above all data points.
#' If cv >= 0, it is also subject to the monotonicity constraint \eqn{\alpha_1 + 2 \alpha_2 t_{j-1} + 3\alpha_3 t_{j-1}^2 + \sum_{l=1}^{j-1} 3 \alpha_{l+3}(t_{j-1}-t_l)^2}.
#' If cv = 1, the concavity constraint \eqn{\pi''(t_j)^T \alpha \leq 0 \quad \forall j = 0,1,...,k_n} is also imposed.
#' The implementation of the cubic spline fitting follows the approach outlined in Daouia, Noh and Park 2016 and uses the
#' framework from the CVXR package.
#' @export
#'
#' @examples
#' xtab=atp_2017$ranking_points
#' ytab=atp_2017$wins
#' xgrid=seq(min(xtab),max(xtab),length.out=101)
#' fit_m <- SplineCubic(xtab=xtab, ytab=ytab, xgrid=xgrid,
#'             kn=2, cv=0, iters=10000, eps=1e-6)
#' fit_mc <- SplineCubic(xtab=xtab, ytab=ytab, xgrid=xgrid,
#'             kn=2, cv=1, iters=100000, eps=1e-6)
#' plot(xtab,ytab, xlab="ranking points", ylab="wins", main="ATP results in 2017 per player")
#' lines(xgrid,fit_m,col="blue")
#' lines(xgrid,fit_mc,col="green")
#' legend("bottomright", col = c("green", "blue"),
#'       legend = c("monontone+concave", "monotone"), lwd = 2)
SplineCubic <- function(xtab, ytab, xgrid, kn, cv, iters=10000, eps=1e-6) {
  x_min <- min(xtab)
  x_max <- max(xtab)
  y_min <- min(ytab)
  y_max <- max(ytab)

  # normalize between 0 and 1
  xtab <- (xtab - x_min) / (x_max - x_min)
  xgrid <- (xgrid - x_min) / (x_max - x_min)
  ytab <- (ytab - y_min) / (y_max - y_min)

  l_end <- min(xtab)
  r_end <- max(xtab)

  # Select FDH or DEA indices (used as knots)
  if (cv > 0) {
    idx <- dea_index(xtab, ytab)
  } else {
    idx <- fdh_index(xtab, ytab)
  }

  ## Interior knots
  t_kn_i <- unname(quantile(xtab[idx], probs = (1:kn) / (kn + 1), type=7)) # type=5 makes it behave like Matlab
  t_kn <- c(l_end, t_kn_i, r_end)

  ## Identity matrix for knots
  n <- length(xtab)
  Ikn <- diag(rep(1, kn + 1))

  ## Build c_j array
  c_j_array <- matrix(0, nrow = kn + 1, ncol = 2 * kn + 3 + 2)

  for (j in 1:(kn + 1)) {
    t_jm1 <- t_kn[j]
    t_j <- t_kn[j + 1]

    if (j == 1) {
      d_j <- c(0, 1, 2 * t_jm1, 6 * t_jm1^2 - 6 * t_jm1 * t_j + 3 * t_j^2, rep(0, kn))
    } else {
      d_j <- c(0, 1, 2 * t_jm1, 6 * t_jm1^2 - 6 * t_jm1 * t_j + 3 * t_j^2)
      for (k in 1:(j - 1)) {
        t_k <- t_kn[k + 1]
        d_j <- c(d_j, 3 * ((t_jm1 - t_k)^2 + (t_j - t_jm1)^2))
      }
      d_j <- c(d_j, rep(0, kn - (j - 1)))
    }
    c_j_array[j, ] <- c(d_j, Ikn[j, ])
  }

  ## Build A_j array
  A_j_array <- array(0, dim = c(2, 2 * kn + 3 + 2, kn + 1))

  for (j in 1:(kn + 1)) {
    t_jm1 <- t_kn[j]
    t_j <- t_kn[j + 1]

    if (j == 1) {
      B_1j <- c(0, 1, 2 * t_jm1, 6 * t_jm1 * t_j - 3 * t_j^2, rep(0, kn))
      B_2j <- c(0, 0, 2 * (t_j - t_jm1), 6 * (t_j - t_jm1) * t_jm1, rep(0, kn))
    } else {
      B_1j <- c(0, 1, 2 * t_jm1, 6 * t_jm1 * t_j - 3 * t_j^2)
      for (k in 1:(j - 1)) {
        t_k <- t_kn[k + 1]
        B_1j <- c(B_1j, 3 * ((t_jm1 - t_k)^2 - (t_j - t_jm1)^2))
      }
      B_1j <- c(B_1j, rep(0, kn - (j - 1)))

      B_2j <- c(0, 0, 2 * (t_j - t_jm1), 6 * (t_j - t_jm1) * t_jm1)
      for (k in 1:(j - 1)) {
        t_k <- t_kn[k + 1]
        B_2j <- c(B_2j, 6 * (t_j - t_jm1) * (t_jm1 - t_k))
      }
      B_2j <- c(B_2j, rep(0, kn - (j - 1)))
    }

    A_j_array[1, , j] <- c(B_1j, -Ikn[j, ])
    A_j_array[2, , j] <- c(B_2j, -Ikn[j, ])
  }

  ## Concavity constraints
  con_array <- do.call(rbind, lapply(0:(kn+1), function(j) {
    x <- t_kn[j + 1]
    c(0, 0, 2, 6 * x, pmax(x - t_kn_i, 0) * 6, rep(0, kn + 1))
  }))

  ## Envelopment constraints
  env_array <- do.call(rbind, lapply(1:n, function(i) {
    x <- xtab[i]
    c(1, x, x^2, x^3, pmax(x - t_kn_i, 0)^3, rep(0, kn + 1))
  }))

  ## Optimization coefficients
  coef <- c(
    (r_end - l_end),
    0.5 * (r_end^2 - l_end^2),
    (1/3) * (r_end^3 - l_end^3),
    (1/4) * (r_end^4 - l_end^4),
    (1/4) * (r_end - t_kn_i)^4,
    rep(0, kn + 1)
  )

  ## Positive parameters constraint
  par_pos <- cbind(matrix(0, nrow = kn + 1, ncol = kn + 3 + 1), diag(rep(1, kn + 1)))

  ## Define optimization variable
  u <- CVXR::Variable(2 * kn + 3 + 2)

  ## Build constraints
  constraints <- list()

  if (cv >= 0) {
  for (j in 1:(kn + 1)) { # monotonicity constraint
    constraints <- c(constraints, CVXR::norm(A_j_array[, , j] %*% u, "2") <= t(c_j_array[j, ]) %*% u)
  }

  }
  if (cv == 1) {
    constraints <- c(constraints, con_array %*% u <= 0)
  }

  constraints <- c(constraints,
                   env_array %*% u >= ytab,
                   par_pos %*% u >= 0,
                   t(coef) %*% u >= 0)

  ## Solve optimization
  objective <- CVXR::Minimize(t(coef) %*% u)
  problem <- CVXR::Problem(objective, constraints)
  result <- CVXR::solve(problem, solver = "SCS", eps = eps, max_iters = iters)

  if (result$status != "optimal") {
    print(paste("DCP problem not optimally solved with k =",k," and iters =",iters,", maybe increase iterations or decrease k"))
  }
  else {
  u_opt <- result$getValue(u)

  alpha <- u_opt[1:(kn + 3 + 1)]

  ## Evaluate fitted spline
  r <- length(xgrid)
  plot_array <- matrix(0, nrow = r, ncol = kn + 3 + 1)

  for (i in 1:r) {
    xpt <- xgrid[i]
    plot_array[i, ] <- c(1, xpt, xpt^2, xpt^3, pmax(xpt - t_kn_i, 0)^3)
  }

  fitt <- plot_array %*% alpha
  fitt <- fitt * (y_max - y_min) + y_min

  return(fitt)
  }
}


#' Fit Upper Frontier Cubic Spline with Optimal Number of Knots
#'
#' fits an (un)constrained cubic spline as upper frontier, whose number of knots is chosen optimally by an information criterion.
#'
#' @param xtab numeric vector of x-indices of data points
#' @param ytab numeric vector of y-indices of data points
#' @param xgrid numeric vector of x-values where function should be evaluated
#' @param cv integer. indicates constraints for fitting: -1 (unconstrained), 0 (monotonically increasing), 1 (monotonically increasing + concave)
#' @param ic string. information criterion: either "AIC" or "BIC"
#' @param kmax integer. highest number of knots considered by information criterion
#' @param iters integer. iterations used for optimization for the CVXR solver
#' @param eps integer. stopping criterion (absolute some of deviations from constraints)
#' @return A list with:
#' \describe{
#'   \item{\code{fit}}{Numeric vector of fitted y-values at the grid points.}
#'   \item{\code{k}}{The selected number of knots.}
#' }
#' @details This function fits cubic splines with different numbers of knots (from 1 to the kmax chosen)
#' and then selects the one with the lowest AIC or BIC value. The information criteria are calculated as
#' \deqn{AIC = \log(\sum|y-i-\hat{phi}(x_i)|) + 2 \cdot (k+4) / n} and
#' \deqn{BIC = \log(\sum|y-i-\hat{phi}(x_i)|) + \log(n) \cdot (k+4) / n}
#' @export
#'
#' @examples
#' xtab = atp_2017$ranking_points
#' ytab = atp_2017$wins
#' xgrid = seq(min(xtab),max(xtab),length.out=101)
#' fit_m <- SplineCubicAuto(xtab=xtab, ytab=ytab, xgrid=xgrid,
#'                 cv=0, ic="BIC", kmax=3, iters=1000000, eps=1e-6)
#' fit_mc <- SplineCubicAuto(xtab=xtab, ytab=ytab, xgrid=xgrid,
#'                 cv=1, ic="BIC", kmax=3, iters=1000000, eps=1e-6)
#' plot(xtab,ytab, xlab="ranking points", ylab="wins", main="ATP results in 2017 per player")
#' lines(xgrid,fit_m$fit,col="red")
#' lines(xgrid,fit_mc$fit,col="blue")
#' legend("bottomright", col = c("blue", "red"),
#'         legend = c("monontone+concave", "monotone"), lwd = 2)
SplineCubicAuto <- function(xtab, ytab, xgrid, cv, ic="BIC", kmax=5, iters=10000, eps=1e-6) {

  if (ic == "AIC") {
    AIC <- numeric()  # Initialize an empty vector for AIC values

    # Loop over k from 1 to 10
    for (k in 1:kmax) {
      # Calculate AIC for current k
      fit <- SplineCubic(xtab, ytab, xtab, k, cv, iters, eps)  # Call SplineCubic function
      residuals <- abs(ytab - fit)  # Calculate residuals
      AIC[k] <- log(sum(residuals)) + 2 * (k + 4) / length(xtab)  # Store AIC value
    }

    # Find the k value that minimizes AIC
    k <- which.min(AIC)  # Find the index of the minimum AIC
    fit <- SplineCubic(xtab, ytab, xgrid, k, cv, iters, eps)  # Call SplineCubic with optimal k

    # Return the fitted values and optimal k
    return(list(fit = fit, k = k))
  }
  if (ic == "BIC") {
    BIC <- numeric()  # Initialize an empty vector for BIC values

    # Loop over k from 1 to 10
    for (k in 1:kmax) {
      # Calculate BIC for current k
      fit <- SplineCubic(xtab, ytab, xtab, k, cv, iters, eps)  # Call SplineCubic function
      residuals <- abs(ytab - fit)  # Calculate residuals
      BIC[k] <- log(sum(residuals)) + log(length(xtab)) * (k + 4) / (length(xtab))  # Store BIC value
    }

    # Find the k value that minimizes BIC
    k <- which.min(BIC)  # Find the index of the minimum BIC
    fit <- SplineCubic(xtab, ytab, xgrid, k, cv, iters, eps)  # Call SplineCubic with optimal k

    # Return the fitted values and optimal k
    return(list(fit = fit, k = k))
  }
}


#' Fit Lower Frontier Cubic Spline with Optimal Number of Knots
#'
#' fits an (un)constrained cubic spline as lower frontier, whose number of knots is chosen optimally by an information criterion.
#'
#' @param xtab numeric vector of x-indices of data points
#' @param ytab numeric vector of y-indices of data points
#' @param xgrid numeric vector of x-values where function should be evaluated
#' @param cv integer. indicates constraints for fitting: -1 (unconstrained), 0 (monotonically increasing), 1 (monotonically increasing + concave)
#' @param ic string. information criterion: either "AIC" or "BIC"
#' @param kmax integer. highest number of knots considered by information criterion
#' @param iters integer. iterations used for optimization for the CVXR solver
#' @param eps integer. stopping criterion (absolute some of deviations from constraints)
#' @return A list with:
#' \describe{
#'   \item{\code{fit}}{Numeric vector of fitted y-values at the grid points.}
#'   \item{\code{k}}{The selected number of knots.}
#' }
#' @details This function fits lower cubic splines with different numbers of knots (from 1 to the kmax chosen)
#' and then selects the one with the lowest AIC or BIC value. The information criteria are calculated as
#' \deqn{AIC = \log(\sum|y-i-\hat{phi}(x_i)|) + 2 \cdot (k+4) / n} and
#' \deqn{BIC = \log(\sum|y-i-\hat{phi}(x_i)|) + \log(n) \cdot (k+4) / n}
#' @export
#'
#' @examples
#' xtab = movies$Year
#' ytab = log(movies$Budget)
#' xgrid = seq(min(xtab),max(xtab),length.out=101)
#' fit_m <- SplineCubicAutoLower(xtab=xtab, ytab=ytab, xgrid=xgrid,
#'                              kmax=5, cv=0, iters=100000, eps=1e-6)
#' plot(xtab,ytab, xlab="year", ylab="log(budget)", main="Runtime of big franchise movies")
#' lines(xgrid,fit_m$fit,col="red")
SplineCubicAutoLower <- function(xtab, ytab, xgrid, cv, ic="BIC", iters=10000, kmax=5, eps = 1e-6) {

  if (ic=="AIC") {
    AIC <- numeric()

    for (k in 1:kmax) {
      fit <- SplineCubicLower(xtab, ytab, xtab, k, cv, iters)
      resid <- ytab - fit
      AIC_k <- log(sum(abs(resid))) + 2 * (k + 4) / length(xtab)
      AIC <- c(AIC, AIC_k)
    }

    k <- which.min(AIC)
    fit <- SplineCubicLower(xtab, ytab, xgrid, k, cv, iters, eps)

    return(list(fit = fit, k = k))
  }

  if (ic == "BIC") {
    BIC <- numeric()  # Initialize an empty vector for BIC values

    for (k in 1:kmax) {
      fit <- SplineCubicLower(xtab, ytab, xtab, k, cv, iters, eps)
      resid <- ytab - fit
      BIC_k <- log(sum(abs(resid))) + log(length(xtab)) * (k + 4) / (2*length(xtab))
      BIC <- c(BIC, BIC_k)
    }

    k <- which.min(BIC)
    fit <- SplineCubicLower(xtab, ytab, xgrid, k, cv, iters, eps)

    return(list(fit = fit, k = k))
  }
}




#' Fit Lower Frontier Cubic Spline with Predefined Number of Knots
#'
#' fits an (un)constrained cubic spline as lower frontier with predefined number of knots
#'
#' @param xtab numeric vector of x-indices of data points
#' @param ytab numeric vector of y-indices of data points
#' @param xgrid numeric vector of x-values where function should be evaluated
#' @param kn integer. number of knots for spline fitting
#' @param cv integer. indicates constraints for fitting: -1 (unconstrained), 0 (monotonically increasing), 1 (monotonically increasing + concave)
#' @param iters integer. iterations used for optimization for the CVXR solver
#' @param eps integer. stopping criterion (absolute some of deviations from constraints)
#' @return numeric vector with y-values of fitted cubic spline at grid points
#' @details The lower boundary curve fitted by this function is a cubic spline with selected number of knots \eqn{t_j} (location chosen at quantiles of $x$)
#' and the solution to the optimization problem \deqn{\min_\alpha \int \pi(x)^T \alpha dx = \sum_{j=0}^{N-1} \alpha_j \int \pi_j(x) dx}
#' subject to the envelopment constraint \eqn{\pi(x_i)^T \alpha \leq y_i \forall i=1,...,n}.
#' Intuitively, the boundary is chosen to maximize the area under the fitted curve while also being weakly below all data points.
#' If cv >= 0, it is also subject to the monotonicity constraint \eqn{\alpha_1 + 2 \alpha_2 t_{j-1} + 3\alpha_3 t_{j-1}^2 + \sum_{l=1}^{j-1} 3 \alpha_{l+3}(t_{j-1}-t_l)^2}.
#' If cv = 1, the concavity constraint \eqn{\pi''(t_j)^T \alpha \leq 0 \quad \forall j = 0,1,...,k_n} is also imposed.
#' The implementation of the cubic spline fitting follows the approach outlined in Daouia, Noh and Park 2016 and uses the
#' framework from the CVXR package.
#' @export
#'
#' @examples
#' xtab=movies$Year
#' ytab=log(movies$Budget)
#' xgrid=seq(min(xtab),max(xtab),length.out=101)
#' fit_m <- SplineCubicLower(xtab=xtab, ytab=ytab, xgrid=xgrid,
#'             kn=2, cv=0, iters=100000, eps=1e-6)
#' fit_mc <- SplineCubicLower(xtab=xtab, ytab=ytab, xgrid=xgrid,
#'             kn=2, cv=1, iters=100000, eps=1e-6)
#' plot(xtab,ytab, xlab="year", ylab="log(budget)", main="Budget of big franchise movies")
#' lines(xgrid,fit_m,col="blue")
#' lines(xgrid,fit_mc,col="green")
#' legend("topright", col = c("green", "blue"),
#'       legend = c("monontone+concave", "monotone"), lwd = 2)
SplineCubicLower <- function(xtab, ytab, xgrid, kn, cv, iters = 10000, eps = 1e-6) {
  x_min <- min(xtab)
  x_max <- max(xtab)
  y_min <- min(ytab)
  y_max <- max(ytab)

  # normalize between 0 and 1
  xtab <- (xtab - x_min) / (x_max - x_min)
  xgrid <- (xgrid - x_min) / (x_max - x_min)
  ytab <- (ytab - y_min) / (y_max - y_min)

  l_end <- min(xtab)
  r_end <- max(xtab)

  if (cv > 0) {
      idx <- dea_index(xtab, ytab)
  } else {
      idx <- fdh_index(xtab, ytab)
  }

  n <- length(xtab)

  # Interior knots
  t_kn_i <- quantile(xtab[idx], probs = (1:kn) / (kn + 1), type=7)
  t_kn <- c(l_end, t_kn_i, r_end)

  Ikn <- diag(kn + 1)
  c_j_array <- list()

  for (j in 1:(kn + 1)) {
    t_jm1 <- t_kn[j]
    t_j <- t_kn[j + 1]
    if (j == 1) {
      d_j <- c(0, 1, 2 * t_jm1, 6 * t_jm1^2 - 6 * t_jm1 * t_j + 3 * t_j^2, rep(0, kn))
    } else {
      d_j <- c(0, 1, 2 * t_jm1, 6 * t_jm1 * t_j - 3 * t_j^2)
      for (k in 1:(j - 1)) {
        t_k <- t_kn[k + 1]
        d_j <- c(d_j, 3 * ((t_jm1 - t_k)^2 - (t_j - t_jm1)^2))
      }
      d_j <- c(d_j, rep(0, kn - (j - 1)))
    }
    c_j <- c(d_j, Ikn[j, ])
    c_j_array[[j]] <- c_j
  }

  A_j_array <- array(0, dim = c(2, 2 * kn + 3 + 2, kn + 1))

  for (j in 1:(kn + 1)) {
    t_jm1 <- t_kn[j]
    t_j <- t_kn[j + 1]
    if (j == 1) {
      B_1j <- c(0, 1, 2 * t_jm1, 6 * t_jm1 * t_j - 3 * t_j^2, rep(0, kn))
      B_2j <- c(0, 0, 2 * (t_j - t_jm1), 6 * (t_j - t_jm1) * t_jm1, rep(0, kn))
    } else {
      B_1j <- c(0, 1, 2 * t_jm1, 6 * t_jm1 * t_j - 3 * t_j^2)
      for (k in 1:(j - 1)) {
        t_k <- t_kn[k + 1]
        B_1j <- c(B_1j, 3 * ((t_jm1 - t_k)^2 - (t_j - t_jm1)^2))
      }
      B_1j <- c(B_1j, rep(0, kn - (j - 1)))

      B_2j <- c(0, 0, 2 * (t_j - t_jm1), 6 * (t_j - t_jm1) * t_jm1)
      for (k in 1:(j - 1)) {
        t_k <- t_kn[k + 1]
        B_2j <- c(B_2j, 6 * (t_j - t_jm1) * (t_jm1 - t_k))
      }
      B_2j <- c(B_2j, rep(0, kn - (j - 1)))
    }
    A_1j <- c(B_1j, -Ikn[j, ])
    A_2j <- c(B_2j, -Ikn[j, ])
    A_j_array[,,j] <- rbind(A_1j, A_2j)
  }

  con_array <- matrix(0, nrow = kn + 2, ncol = 2 * kn + 3 + 2)
  for (j in 0:(kn + 1)) {
    x_val <- t_kn[j + 1]
    pitwo <- c(0, 0, 2, 6 * x_val, pmax(0, x_val - t_kn_i) * 6)
    con_array[j + 1, ] <- c(pitwo, rep(0, kn + 1))
  }

  env_array <- matrix(0, nrow = n, ncol = 2 * kn + 3 + 2)
  for (i in 1:n) {
    x_val <- xtab[i]
    pizero <- c(1, x_val, x_val^2, x_val^3, pmax(0, x_val - t_kn_i)^3)
    env_array[i, ] <- c(pizero, rep(0, kn + 1))
  }

  coef <- c(
    r_end - l_end,
    0.5 * (r_end^2 - l_end^2),
    (1/3) * (r_end^3 - l_end^3),
    (1/4) * (r_end^4 - l_end^4),
    (1/4) * (r_end - t_kn_i)^4,
    rep(0, kn + 1)
  )

  u <- CVXR::Variable(2 * kn + 3 + 2)
  constraints <- list()

  if (cv >= 0) {
  for (j in 1:(kn + 1)) {
    A_j <- A_j_array[,,j]
    c_j <- c_j_array[[j]]
    constraints <- c(constraints, CVXR::norm(A_j %*% u, "2") <= sum(c_j * u))
      # monotonicity constraint
  }
  }
  if (cv == 1) {
    constraints <- c(constraints, con_array %*% u >= 0)
  }

  constraints <- c(
    constraints,
    env_array %*% u <= ytab,
    diag(kn + 1) %*% u[(2 * kn + 3 + 2 - kn):(2 * kn + 3 + 2)] >= 0,
    sum(coef * u) >= 0
  )

  objective <- CVXR::Maximize(sum(coef * u))
  problem <- CVXR::Problem(objective, constraints)
  result <- CVXR::solve(problem, solver = "SCS", eps = eps, max_iters = iters)

  if (result$status == "optimal") {
    alpha <- result$getValue(u)[1:(kn + 3 + 1)]

    plot_array <- t(sapply(xgrid, function(xpt) {
      c(1, xpt, xpt^2, xpt^3, pmax(0, xpt - t_kn_i)^3)
    }))

    fitt <- plot_array %*% alpha
    fitt <- fitt * (y_max - y_min) + y_min

    return(fitt)
  } else {
    print(paste("DCP problem not optimally solved with k =",k,", eps = ",eps," and iters =",iters,", maybe increase iterations or decrease k"))
  }
}
