/*
    Outcome Model from 'cookstove' R package

Model Notes:
    --Exposure coefficients *vary* across all studies
    --Beta restriction
    --No cluster-level RE
*/
data {
    // Dimensions
    int<lower=1> S; //total number of studies
    int<lower=1> N; //total number of observations
    int<lower=1> n; //total number of subjects
    int<lower=1> xdf; // number of basis functions for exposure (e.g. PM)
    int<lower=0> timedf; // degrees of freedom in time spline
    int<lower=0> p; // number of covariates -- do NOT include intercept here
    // Data
    matrix[p==0 ? 0 : N,p] Z; // Matrix of covariates (not counting exposure, time, or intercept)
    matrix[N,xdf] Hx; // Matrix of exposure values, transformed to basis
    matrix[timedf==0 ? 0 : N,timedf] Ht;  // matrix of time spline values
    int<lower=0> y[N]; // Case indicator
    int<lower=0> nT[N]; // Time at risk
    int<lower=1, upper=n> subj_of_obs[N]; // which subj is this obs?
    int<lower=1, upper=S> study_of_subj[n]; // which study is this subj in?
    // Hyperparameters
    real prior_sigmaI_mean;  // Prior for SD of subj-level RE
    real<lower=0> prior_sigmaI_sd;
    real prior_delta_mean; // Prior mean for coefs delta (usually zero)
    real prior_sigmaDelta_mean; // Prior for SD of coefs delta
    real<lower=0> prior_sigmaDelta_sd;
    real prior_gamma_mean; // Prior mean for coefs gamma (usually zero)
    real prior_sigmaGamma_mean; // Prior for SD of coefs gamma
    real<lower=0> prior_sigmaGamma_sd;
    real prior_beta_mean; // Prior mean for coefs beta (usually zero)
    real prior_sigmaBeta_mean; // Prior for SD of coefs beta
    real<lower=0> prior_sigmaBeta_sd;
    real<lower=0> beta_nu;
}
transformed data {
    int<lower=1, upper=S> study_of_obs[N]; // which study is this subj in?
    study_of_obs=study_of_subj[subj_of_obs];
}
parameters {
       real<lower=0> sigmaI;
    real<lower=0> sigmaGamma[p==0 ? 0 : 1] ;
    real<lower=0> sigmaDelta[timedf==0 ? 0 : 1] ;
    real<lower=0> sigmaBeta;
    vector[S] bS; // study-level intercept
    vector[p] gamma_raw; // coefficients for covariates Z
    vector[timedf] delta_raw; // coefficients for time spline Ht
    matrix[S, xdf] beta_raw; // Coefficients for exposure spline
    vector[n] reI_raw;
    cholesky_factor_corr[S] L;
}
transformed parameters {
    vector[N] mui; // Mean of observations
    vector[n] reI; // random effects by subject
    vector[p] gamma; // coefficients for covariates Z
    vector[timedf] delta; // coefficients for time spline Ht
    matrix[S, xdf] beta; // Coefficients for exposure spline
    vector[N] Hxbeta;
    reI = sigmaI*reI_raw;
    beta = prior_beta_mean + L*(sigmaBeta *beta_raw);
    Hxbeta = rows_dot_product(Hx, beta[study_of_obs]);
    mui= bS[study_of_obs] + Hxbeta + reI[subj_of_obs];
    if (p >0 ){
        gamma = prior_gamma_mean + sigmaGamma[1]*gamma_raw;
        mui += Z*gamma;
    }
    if (timedf > 0){
        delta = prior_delta_mean + sigmaDelta[1]*delta_raw;
        mui += Ht*delta;
    }
}
model {
    if (p >0) {
        target += normal_lpdf(gamma_raw | 0, 1);
        target += normal_lpdf(sigmaGamma | prior_sigmaGamma_mean, prior_sigmaGamma_sd);
    }
    if (timedf > 0){
        target += normal_lpdf(sigmaDelta | prior_sigmaDelta_mean, prior_sigmaDelta_sd);
        target += normal_lpdf(delta_raw | 0, 1);
    }
    target += normal_lpdf(to_vector(beta_raw) | 0, 1);
    target += lkj_corr_cholesky_lpdf(L | beta_nu);
    target += normal_lpdf(reI_raw | 0, 1);
    target += normal_lpdf(sigmaI | prior_sigmaI_mean, prior_sigmaI_sd);
    target += normal_lpdf(sigmaBeta | prior_sigmaBeta_mean, prior_sigmaBeta_sd);
    target += binomial_logit_lpmf(y | nT, mui);
}
generated quantities {
    matrix[S,S] betaCorr;
    betaCorr = multiply_lower_tri_self_transpose(L);
}


