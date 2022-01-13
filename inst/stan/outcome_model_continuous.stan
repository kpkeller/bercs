/*
    Outcome Model for 'bercs' R package
*/
data {
    // Dimensions
// Outcome model dimensions

int<lower=1> S; //total number of studies
int<lower=1> N; //total number of observations
int<lower=1> n; //total number of units
int<lower=0> xdf; // number of basis functions for exposure (e.g. PM)
int<lower=0> timedf; // degrees of freedom in time spline
int<lower=0> p; // number of covariates -- do NOT include intercept here
    // Data
// Data statements for outcome model
matrix[p==0 ? 0 : N,p] Z; // Matrix of covariates (not counting exposure, time, or intercept)
matrix[N,xdf] Mx; // Matrix of exposure values, transformed to basis
matrix[timedf==0 ? 0 : N,timedf] Mt;  // matrix of time spline values
real y[N]; // Case indicator
int<lower=0> nT[N]; // Time at risk
int<lower=1, upper=n> unit_of_obs[N]; // which unit is this obs?
int<lower=1, upper=S> study_of_unit[n]; // which study is this unit in?

    // declares hyperpriors
// Outcome model hyperparameters
// for reI, delta, gamma

real prior_sigmaI_mean;  // Prior for SD of unit-level RE
real<lower=0> prior_sigmaI_sd;
real prior_delta_mean; // Prior mean for coefs delta (usually zero)
real prior_sigmaDelta_mean; // Prior for SD of coefs delta
real<lower=0> prior_sigmaDelta_sd;
real prior_gamma_mean; // Prior mean for coefs gamma (usually zero)
real prior_sigmaGamma_mean; // Prior for SD of coefs gamma
real<lower=0> prior_sigmaGamma_sd;
//new parameter, sigma_y
real<lower=0> prior_sigma_y_mean; //prior for new continuous sigma mean
real<lower=0> prior_sigma_y_sd; //prior for new continuous sigma sd

    real prior_beta_mean; // Prior mean for coefs beta (usually zero)
    real prior_sigmaBeta_mean; // Prior for SD of coefs beta
    real<lower=0> prior_sigmaBeta_sd;
}
transformed data {
#include /data/outcome_tdata.stan
}
parameters {

real<lower=0> sigmaI;
real<lower=0> sigmaGamma[p==0 ? 0 : 1] ;
real<lower=0> sigmaDelta[timedf==0 ? 0 : 1] ;
vector[S] bS; // Study-level intercept
vector[p] gamma_raw; // Coefficients for covariates Z
vector[timedf] delta_raw; // Coefficients for time spline Ht
vector[n] reI_raw;

    real<lower=0> sigmaBeta;
    vector[xdf] beta_raw; // Coefficients for exposure spline
    //new parameter sigma
    real<lower=0> sigma_y;
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

// Model statement for outcome model
// except for updates to beta and related parameters

if (p >0) {
    target += normal_lpdf(gamma_raw | 0, 1);
    target += normal_lpdf(sigmaGamma | prior_sigmaGamma_mean, prior_sigmaGamma_sd);
}
if (timedf > 0){
    target += normal_lpdf(sigmaDelta | prior_sigmaDelta_mean, prior_sigmaDelta_sd);
    target += normal_lpdf(delta_raw | 0, 1);
}
target += normal_lpdf(reI_raw | 0, 1);
target += normal_lpdf(sigmaI | prior_sigmaI_mean, prior_sigmaI_sd);
target += normal_lpdf(sigma_y | prior_sigma_y_mean, prior_sigma_y_sd);
//added sigma_y to model
target += normal_lpdf(y | mui, sigma_y);

    if (xdf> 0){
        target += normal_lpdf(beta_raw | 0, 1);
        target += normal_lpdf(sigmaBeta | prior_sigmaBeta_mean, prior_sigmaBeta_sd);
    }
}

