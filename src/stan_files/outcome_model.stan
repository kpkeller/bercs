/*
    Outcome Model for 'bercs' R package
*/
data {
    // Dimensions
#include /data/outcome_dim.stan
    // Data
#include /data/outcome_data.stan
    // declares hyperpriors
#include /data/outcome_hyperprior.stan
    real prior_beta_mean; // Prior mean for coefs beta (usually zero)
    real prior_sigmaBeta_mean; // Prior for SD of coefs beta
    real<lower=0> prior_sigmaBeta_sd;
}
transformed data {
#include /data/outcome_tdata.stan
}
parameters {
#include /parameters/parameters_outcome.stan
    real<lower=0> sigmaBeta;
    vector[xdf] beta_raw; // Coefficients for exposure spline
}
transformed parameters {
#include tparameters/tparameters_outcome_define.stan
    vector[xdf] beta; // Coefficients for exposure spline
#include tparameters/tparameters_outcome_update.stan
    if (xdf > 0){
        beta = prior_beta_mean + sigmaBeta*beta_raw;
        mui += Mx*beta;
    }

}
model {
#include /model/model_stem_outcome.stan
    if (xdf> 0){
        target += normal_lpdf(beta_raw | 0, 1);
        target += normal_lpdf(sigmaBeta | prior_sigmaBeta_mean, prior_sigmaBeta_sd);
    }

}



