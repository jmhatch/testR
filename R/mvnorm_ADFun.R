#' Create a \code{TMB::ADFun} object for the multivariate normal likelihood.
#'
#' @param X Matrix of observations.
#' @return A list as returned by \code{TMB::MakeADFun} representing the negative loglikelihood of a multivariate normal.
#' @export
mvnorm_ADFun <- function(X) {
  TMB::MakeADFun(data = list(model = "MVNormalNLL", X = X),
                 parameters = list(mu = rep(0, 2), sd = rep(1, 2), rho = 0.5),
                 DLL = "testR_TMBExports", silent = TRUE)
}
