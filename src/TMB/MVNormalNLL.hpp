/// @file MVNormalNLL.hpp

#ifndef MVNormalNLL_hpp
#define MVNormalNLL_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// Negative log-likelihood of the multivariate normal distribution.
template<class Type>
Type MVNormalNLL(objective_function<Type>* obj) {
  using namespace density;
  matrix<Type> Sigma(3,3);
  Sigma.fill(0.1);             // Fill the whole matrix
  Sigma.diagonal() *= 10.0;    // Multiply diagonal by 10 to positive definite Sigma
  vector<Type> x0(3);          // Point of evaluation
  x0.fill(0.0);                // Initialize x0 to be zero
  MVNORM_t<Type> N_0_Sigma(Sigma);   // N_0_Sigma is now a Distribution
  vector<Type> res(3);
  res = N_0_Sigma(x0);         // Evaluates (neg. log) density at x
  return res.sum();
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
