#--- functions for the generalized hyperbolic regression model -----------------

# model is:
# y_i = x_i'beta + v_i * alpha + exp(w_i'gamma) v_i^(1/2) z_i,
# z_i ~iid N(0,1)
# v_i ~iid GIG(psi, 1, lambda)

#--- generalized inverse gaussian distribution ---------------------------------

# C++ code for rgig
# requires installation of Rcpp package
Rcpp::sourceCpp("randGIG.cpp")

#' Density of the Generalized Inverse-Gaussian (GIG) distribution
#'
#' @param x Vector of quantiles.
#' @param psi Vector of shape parameters.
#' @param eta Vector of scale parameters.
#' @param lambda Vector of power parameters.
#' @param log Logical; whether to return density on log scale.
#' @return Vector of density evaluations.
#' @details If \code{X ~ GIG(psi, eta, lambda)}, then the PDF of \code{X} is
#' \preformatted{
#' f(x | psi, eta, lambda) = cst * x^(lambda-1) * exp(-psi/2 * (x/eta + eta/x)),
#' }
#' for \code{x > 0}.
dgig <- function(x, psi, eta, lambda, log = FALSE) {
  cst <- -log(2*eta*besselK(psi, nu = lambda, expon.scaled = TRUE))+psi
  xe <- x/eta
  ans <- cst + log(xe) * (lambda-1) - .5 * psi * (xe + 1/xe)
  if(!log) ans <- exp(ans)
  ans
}

#' Generate random draws from the GIG distribution.
#'
#' @param n Number of random samples to generate.
#' @param psi Vector of shape parameters.
#' @param eta Vector of scale parameters.
#' @param lambda Vector of power parameters.
#' @return Vector of random draws.
rgig <- function(n, psi, eta, lambda) {
  # format parameter inputs
  if(length(psi) != 1) psi <- rep(psi, len = n)
  if(length(eta) != 1) eta <- rep(eta, len = n)
  if(length(lambda) != 1) lambda <- rep(lambda, len = n)
  .GenerateGIG(n, psi, eta, lambda) # entry point to C++ code
}

#--- generalized hyperbolic distribution ---------------------------------------

#' Density of the Generalized Hyperbolic (GH) distribution.
#'
#' @param x Vector of quantiles.
#' @param mu Vector of location parameters.
#' @param alpha Vector of skewness parameters.
#' @param sigma Vector of scale parameters.
#' @param psi Vector of shape parameters.
#' @param lambda Vector of power parameters.
#' @param log Logical; whether to return density on log scale.
#' @return Vector of density evaluations.
#' @details If \code{X ~ GH(mu, alpha, sigma, psi, lambda)}, then it is generated via
#' \preformatted{
#' V ~ GIG(psi, 1, lambda)
#' X | V ~ N(mu + a V^(1/2), sigma^2 V).
#' }
#' The PDF of \code{X} is given by
#' \preformatted{
#' Q(x)^(lambda/2-1/4) * besselK(sqrt(Q(x) * R), lambda-1/2) * exp((x-mu)alpha/sigma),
#' }
#' where \code{Q(x) = psi + (x-mu)^2/sigma^2} and \code{R = psi + (alpha/sigma)^2}.
dghyp <- function(x, mu, alpha, sigma, psi, lambda, log = FALSE) {
  ls2pi <- 0.918938533204672741780329736406 # log(2*pi)/2
  z <- (x-mu)/sigma # standardize
  oas <- psi + (alpha/sigma)^2
  oz2 <- psi + z^2
  soaz <- sqrt(oas * oz2)
  # all logs together except power
  ans <- log(besselK(psi, nu = lambda, expon.scaled = TRUE) /
             besselK(soaz, nu = lambda - .5, expon.scaled = TRUE) * sigma)
  # add in power = product-log
  ans <- - ans - soaz + psi + .5 * (lambda -.5) * log(oz2/oas)
  ans <- ans + z*alpha/sigma - ls2pi # log-less part
  if(!log) ans <- exp(ans)
  ans
}

#' Generate random draws from the GH distribution.
#'
#' @param n Number of random samples to generate.
#' @param mu Vector of location parameters.
#' @param alpha Vector of skewness parameters.
#' @param sigma Vector of scale parameters.
#' @param psi Vector of shape parameters.
#' @param lambda Vector of power parameters.
#' @return Vector of random draws.
rghyp <- function(n, mu, alpha, sigma, psi, lambda) {
  V <- rgig(n, psi, 1, lambda) # mixing parameter
  mu + V * alpha + sigma * sqrt(V) * rnorm(n) # result
}
