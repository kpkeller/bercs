/*
    Outcome Model for 'bercs' R package
    Different curve per study
*/
data {
    // Dimensions
#include /data/outcome_dim.stan
    // Data
#include /data/outcome_data.stan
    // Hyperparameters
#include /data/outcome_hyperprior.stan
    real prior_beta_mean; // Prior mean for coefs beta (usually zero)
    real prior_sigmaBeta_mean; // Prior for SD of coefs beta
    real<lower=0> prior_sigmaBeta_sd;
    real<lower=0> beta_nu;
}
transformed data {
#include /data/outcome_tdata.stan
}
parameters {
#include /parameters/parameters_outcome.stan
    real<lower=0> sigmaBeta;
    matrix[S, xdf] beta_raw; // Coefficients for exposure spline
    cholesky_factor_corr[S] L;
}
transformed parameters {
#include tparameters/tparameters_outcome_define.stan
    matrix[S, xdf] beta; // Coefficients for exposure spline
#include tparameters/tparameters_outcome_update.stan
    if (xdf > 0){
        beta = prior_beta_mean + L*(sigmaBeta *beta_raw);
        mui += rows_dot_product(Mx, beta[study_of_obs]);
    }
}
model {
#include /model/model_stem_outcome.stan
    target += normal_lpdf(to_vector(beta_raw) | 0, 1);
    target += lkj_corr_cholesky_lpdf(L | beta_nu);
    target += normal_lpdf(sigmaBeta | prior_sigmaBeta_mean, prior_sigmaBeta_sd);
}
generated quantities {
    matrix[S,S] betaCorr;
    betaCorr = multiply_lower_tri_self_transpose(L);
}


