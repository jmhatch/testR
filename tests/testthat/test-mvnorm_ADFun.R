
context("mvnorm_ADFun")


test_that("mvnorm_ADFun calculates correct negloglik, gradient, and hessian", {
  nll_fun <- function(theta, X) {
    # mean 
    mu = c(theta[1], theta[2])
    # Sigma
    Sigma = matrix(c(theta[3] * theta[3], theta[3] * theta[4] * theta[5],
                     theta[3] * theta[4] * theta[5], theta[4] * theta[4]), 
                   nrow = 2, ncol = 2, byrow = TRUE)
    -sum(emdbook::dmvnorm(X, mu, Sigma, log = TRUE))
  }
  nreps <- 20
  for(ii in 1:nreps) {
    # simulate data/parameters
    n <- sample(10:100, 1)
    mu = rnorm(2)
    sd = rexp(2)
    rho = runif(1)
    Sigma = matrix(c(sd[1] * sd[1], sd[1] * sd[2] * rho,
                     sd[2] * sd[1] * rho, sd[2] * sd[2]), 
                   nrow = 2, ncol = 2, byrow = TRUE)
    X = MASS::mvrnorm(n = n, mu = mu, Sigma) 
    # nll/g/h/ with R + numDeriv
    theta = c(mu, sd, rho)
    ll1 <- nll_fun(theta, X)
    gg1 <- numDeriv::grad(func = nll_fun, x = theta, X = X)
    hh1 <- numDeriv::hessian(func = nll_fun, x = theta, X = X)
    # nll/g/h with TMB
    nll_obj <- mvnorm_ADFun(X)
    ll2 <- nll_obj$fn(theta)
    gg2 <- nll_obj$gr(theta)[1,]
    hh2 <- nll_obj$he(theta)
    # check they are identical
    expect_equal(ll1, ll2)
    expect_equal(gg1, gg2)
    expect_equal(hh1, hh2)
  }
})
