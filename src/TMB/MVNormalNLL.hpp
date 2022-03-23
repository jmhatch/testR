/// @file MVNormalNLL.hpp (from admb-project/tmb-examples)

#ifndef MVNormalNLL_hpp
#define MVNormalNLL_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// Negative log-likelihood of the multivariate normal distribution.
template<class Type>
Type MVNormalNLL(objective_function<Type>* obj) {
  
  // data
  DATA_MATRIX(X);
  int n = X.rows();
  int p = X.cols();
  
  // parameters
  PARAMETER_VECTOR(mu);
  PARAMETER_VECTOR(sd);
  PARAMETER(rho);
  
  // covariance
  matrix<Type> Sigma(p, p);
  Sigma(0, 0) = sd(0) * sd(0);
  Sigma(0, 1) = sd(0) * sd(1) * rho;
  Sigma(1, 1) = sd(1) * sd(1);
  Sigma(1, 0) = sd(0) * sd(1) * rho;
  
  // allocate space
  vector<Type> residual(p);
  
  // accumulator
  Type neglogL = 0.0;
  
  // mvnormal nll
  using namespace density;
  MVNORM_t<Type> neg_log_dmvnorm(Sigma);
  
  // loop thru data
  for(int i = 0; i < n; i++)
  {
    residual = vector<Type>(X.row(i)) - mu;
    neglogL += neg_log_dmvnorm(residual);
  }
  
  REPORT(Sigma);
  REPORT(sd);
  REPORT(rho);
  
  return neglogL;
  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
