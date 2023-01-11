#' Get vector average by a variable
#'
#' @param vec A vector of numerics
#' @param var A grouping variable
#' @param mode Either 'collapse' or 'sub-sample' 
#' @param n A vector of numerics
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c('group1', 'group2'), 3)
#' n <- sample(2:3, 6, replace = TRUE)
#' pcr:::.pcr_average(vec, var)
#' pcr:::.pcr_average(vec, var, mode = 'sub-sample', n = n)
#'
#' @importFrom stats aggregate

.pcr_average <- function(vec, var, mode = 'collapse', n = NULL) {
  if(mode == 'collapse'){
    res <- aggregate(vec,
                     by = list(var),
                     FUN = mean)
    res <- res[order(unique(var)),]
    return(res$x)
  }
  else if(mode == 'sub-sample'){
    res <- c()
    for (group in unique(var)){
      i <- which(var %in% group)
      mu_aggr <- sum(n[i]*vec[i])/sum(n[i])
      res <- c(res, mu_aggr)
    }
    res <- res[order(unique(var))]
    return(res)
    
    
  }
  else {
    stop("mode can be one of 'collapse' or 'sub-sample'.")
  }
}

#' Helper to get vector aggregated, Bessel-corrected standard deviation 
#' by a sub-sampled variable
#' 
#' @param mu A vector of numerics
#' @param n A vector of numerics
#' @param sd A vector of numerics
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' mu <- rnorm(6, 30, 1)
#' n <- sample(2:3, 6, replace = TRUE)
#' sd <- rnorm(6, .1, .01)
#' pcr:::.pcr_sscv(mu, n, sd)
#'
#' @importFrom stats aggregate sd
#' 
.pcr_aggrsd <- function(mu, n, sd){
  #Sample based mean
  mu_aggr <- sum(n*mu)/sum(n)
  #Sample based SD
  sd_aggr <- sqrt((sum((n-1)*sd^2)+sum(n*(mu-mu_aggr)^2))/(sum(n)-1))
  return(sd_aggr)
}

#' Get vector standard deviation by a variable
#'
#' @inheritParams .pcr_average
#' @param sd A vector of numerics
#' 
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c('group1', 'group2'), 3)
#' n <- sample(2:3, 6, replace = TRUE)
#' sd <- rnorm(6, .1, .01)
#' pcr:::.pcr_sd(vec, var)
#' pcr:::.pcr_sd(vec, var, mode = 'sub-sample', n = n, sd = sd)
#'
#' @importFrom stats aggregate sd

.pcr_sd <- function(vec, var) {
if (mode == 'collapse'){
    res <- aggregate(vec,
                     by = list(var),
                     FUN = sd)
    res <- res[order(unique(var)),] 
    return(res$x)
  }
  else if (mode == 'sub-sample'){
    res <- c()
    for (group in unique(var)){
      i <- which(var %in% group)
      sd_aggr <- .pcr_aggrsd(vec[i], n[i], sd[i])
      res <- c(res, sd_aggr)
    }
    res <- res[order(unique(var))]
    return(res)
  }
  else {
    stop("mode can be one of 'collapse' or 'sub-sample'.")
  }
}

#' Get vector coefficient of variance by a variable
#'
#' @inheritParams .pcr_sd
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c('group1', 'group2'), 3)
#' n <- sample(2:3, 6, replace = TRUE)
#' sd <- rnorm(6, .1, .01)
#' pcr:::.pcr_cv(vec, var)
#' pcr:::.pcr_cv(vec, var, mode = 'sub-sample', n = n, sd = sd)
#'
#' @importFrom stats aggregate sd

.pcr_cv <- function(vec, var, mode = 'collapse', n = NULL, sd = NULL) {
  if (mode == 'collapse'){
    res <- aggregate(vec,
                     by = list(var),
                     FUN = function(x) sd(x)/mean(x))
    res <- res[order(unique(var)),] 
    return(res$x)
  }
  else if (mode == 'sub-sample'){
    res <- c()
    for (group in unique(var)){
      i <- which(var %in% group)
      cv_aggr <- .pcr_aggrsd(vec[i], n[i], sd[i])/(sum(n[i]*vec[i])/sum(n[i]))
      res <- c(res, cv_aggr)
    }
    res <- res[order(unique(var))]
    return(res)
  }
  else {
    stop("mode can be one of 'collapse' or 'sub-sample'.")
  }
}

#' Normalize vector by another
#'
#' @inheritParams .pcr_average
#' @param ref A numeric vector
#' @param mode Either 'subtract' or 'divide'
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' ref <- rnorm(6, 30, .1)
#' pcr:::.pcr_normalize(vec, ref)
#' pcr:::.pcr_normalize(vec, ref, mode = 'divide')
#'
#' @importFrom stats aggregate sd

.pcr_normalize <- function(vec, ref, mode = 'subtract') {
  if (mode == 'subtract') {
    res <- vec - ref
  } else if (mode == 'divide') {
    res <- vec / ref
  } else {
    stop("mode can be one of 'subtract' or 'divide'.")
  }
  return(res)
}

#' Propage two vectors
#'
#' @inheritParams .pcr_average
#' @inheritParams .pcr_normalize
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' ref <- rnorm(6, 30, .1)
#' pcr:::.pcr_error(vec, ref)

.pcr_error <- function(vec, ref) {
  res <- sqrt(vec^2 + ref^2)
  return(res)
}

#' Calculate the amounts
#'
#' @inheritParams .pcr_average
#' @param a A numeric
#' @param b A numeric
#'
#' @return A vector of numerics
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' pcr:::.pcr_amount(vec, 1, 1)

.pcr_amount <- function(vec, a, b) {
  res <- 10 ^ ((vec - a)/b)
  return(res)
}

#' Raise two to a vector power
#'
#' @inheritParams .pcr_average
#'
#' @return A numeric
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' pcr:::.pcr_relative(vec)

.pcr_relative <- function(vec) {
  res <- 2 ^ (-vec)
  return(res)
}

#' Calculate R squared
#'
#' @inheritParams .pcr_average
#'
#' @return A numeric
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c(.1, .5), 3)
#' pcr:::.pcr_rsquared(vec, var)
#'
#' @importFrom stats cor

.pcr_rsquared <- function(vec, var) {
  if(anyNA(vec)){
    warning(paste0(sum(is.na(vec)),
                   " NAs detected. ",
                   "Ensure samples are still in the dynamic range"))
  }
  res <- cor(vec,
             log10(var),
             use = "complete.obs")^2
  return(res)
}

#' Calculate the intercept of a line
#'
#' @inheritParams .pcr_average
#'
#' @return A numeric
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c(.1, .5), 3)
#' pcr:::.pcr_intercept(vec, var)
#'
#' @importFrom stats lm coefficients

.pcr_intercept <- function(vec, var) {
  ll <- lm(vec ~ log10(var))
  res <- coefficients(ll)
  return(res[[1]])
}

#' Calculate the slope of a line
#'
#' @inheritParams .pcr_average
#'
#' @return A numeric
#'
#' @keywords internal
#'
#' @examples
#' vec <- rnorm(6, 30, 1)
#' var <- rep(c(.1, .5), 3)
#' pcr:::.pcr_slope(vec, var)
#'
#' @importFrom stats lm coefficients

.pcr_slope <- function(vec, var) {
  ll <- lm(vec ~ log10(var))
  res <- coefficients(ll)
  return(res[[2]])
}
