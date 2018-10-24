/*
    Outcome Model from 'cookstove' R package

Model Notes:
    --Exposure coefficients shared across all studies
    --No cluster-level RE
*/
data {
    // Dimensions
    int<lower=1> S; //total number of studies
    int<lower=1> N; //total number of observations
    int<lower=1> n; //total number of subjects
    int<lower=0> xdf; // number of basis functions for exposure (e.g. PM)
    int<lower=0> timedf; // degrees of freedom in time spline
    int<lower=0> p; // number of covariates -- do NOT include intercept here
    // Data
    matrix[p==0 ? 0 : N,p] Z; // Matrix of covariates (not counting exposure, time, or intercept)
    matrix[xdf==0 ? 0 : N,xdf] Hx; // Matrix of exposure values, transformed to basis
    matrix[timedf==0 ? 0 : N,timedf] Ht;  // matrix of time spline values
    int<lower=0> y[N]; // Case indicator
    int<lower=0> nT[N]; // Time at risk
    int<lower=1, upper=n> subj_of_obs[N]; // which subj is this obs?
    int<lower=1, upper=S> study_of_subj[n]; // which study is this subj in?
    // Hyperparameters
    real prior_sigmaI_mean;  // Prior for SD of subj-level RE
    real<lower=0> prior_sigmaI_sd;
    real prior_delta_mean; // Prior mean for coefs delta (usually zero)
    real prior_sigmaDel_mean; // Prior for SD of coefs delta
    real<lower=0> prior_sigmaDel_sd;
    real prior_gamma_mean; // Prior mean for coefs gamma (usually zero)
    real prior_sigmaGam_mean; // Prior for SD of coefs gamma
    real<lower=0> prior_sigmaGam_sd;
    real prior_beta_mean; // Prior mean for coefs beta (usually zero)
    real prior_sigmaBeta_mean; // Prior for SD of coefs beta
    real<lower=0> prior_sigmaBeta_sd;
}
transformed data {
    int<lower=1, upper=S> study_of_obs[N]; // which study is this subj in?
    study_of_obs=study_of_subj[subj_of_obs];
}
parameters {
    real<lower=0> sigmaI;
    real<lower=0> sigmaGam[p==0 ? 0 : 1] ;
    real<lower=0> sigmaDel[timedf==0 ? 0 : 1] ;
    real<lower=0> sigmaBeta[xdf==0 ? 0 : 1] ;
    vector[S] bS; // Study-level intercept
    vector[p] gamma_raw; // Coefficients for covariates Z
    vector[timedf] delta_raw; // Coefficients for time spline Ht
    vector[xdf] beta_raw; // Coefficients for exposure spline
    vector[n] reI_raw;
}
transformed parameters {
    vector[N] mui; // Mean of observations
    vector[n] reI; // Random effects by subject
    vector[p] gamma; // Coefficients for covariates Z
    vector[timedf] delta; // Coefficients for time spline Ht
    vector[xdf] beta; // Coefficients for exposure spline
    reI = sigmaI*reI_raw;
    mui= bS[study_of_obs] + reI[subj_of_obs];
    if (xdf > 0){
        beta = prior_beta_mean + sigmaBeta[1]*beta_raw;
        mui += Hx*beta;
    }
    if (p >0 ){
        gamma = prior_gamma_mean + sigmaGam[1]*gamma_raw;
        mui += Z*gamma;
    }
    if (timedf > 0){
            delta = prior_delta_mean + sigmaDel[1]*delta_raw;
        mui += Ht*delta;
    }
}
model {
    if (p >0) {
        target += normal_lpdf(gamma_raw | 0, 1);
        target += normal_lpdf(sigmaGam | prior_sigmaGam_mean, prior_sigmaGam_sd);
    }
    if (timedf > 0){
        target += normal_lpdf(sigmaDel | prior_sigmaDel_mean, prior_sigmaDel_sd);
        target += normal_lpdf(delta_raw | 0, 1);
    }
    if (xdf> 0){
        target += normal_lpdf(beta_raw | 0, 1);
        target += normal_lpdf(sigmaBeta | prior_sigmaBeta_mean, prior_sigmaBeta_sd);
    }
    target += normal_lpdf(reI_raw | 0, 1);
    target += normal_lpdf(sigmaI | prior_sigmaI_mean, prior_sigmaI_sd);
    target += binomial_logit_lpmf(y | nT, mui);
}



