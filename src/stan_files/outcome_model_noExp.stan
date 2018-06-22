/*
    Outcome Model from 'cookstove' R package

Model Notes:
    --No Exposure
    --No cluster-level RE
*/
data {
    // Dimensions
    int<lower=1> N; //total number of observations
    int<lower=1> n; //total number of subjects
    int<lower=1> timedf; // degrees of freedom in time spline
    int<lower=1> p; // number of covariates -- do NOT include intercept here
    // Data
    matrix[N,p] Z; // Matrix of covariates (not counting exposure, time, or intercept)
    matrix[N,timedf] Ht;  // matrix of time spline values
    int<lower=0> y[N]; // Case indicator
    int<lower=0> nT[N]; // Time at risk 
    int<lower=1, upper=n> subj_of_obs[N]; // which subj is this obs?
    // Hyperparameters 
    real prior_sigmaI_mean;  // Prior for SD of subj-level RE
    real<lower=0> prior_sigmaI_sd; 
    real prior_delta_mean; // Prior mean for coefs delta (usually zero)
    real prior_sigmaDel_mean; // Prior for SD of coefs delta
    real<lower=0> prior_sigmaDel_sd; 
    real prior_gamma_mean; // Prior mean for coefs gamma (usually zero)
    real prior_sigmaGam_mean; // Prior for SD of coefs gamma
    real<lower=0> prior_sigmaGam_sd; 
}
parameters {
    real bS; // study-level intercept 
    vector[p] gamma_raw; // Coefficients for covariates Z
    vector[timedf] delta_raw; // Coefficients for time spline Ht
    vector[n] reI_raw;
    real<lower=0> sigmaI;
    real<lower=0> sigmaGam;
    real<lower=0> sigmaDel;
}
transformed parameters {
    vector[N] mui; // Mean of observations
    vector[n] reI; // Random effects by subject
    vector[p] gamma; // Coefficients for covariates Z
    vector[timedf] delta; // Coefficients for time spline Ht
    reI = sigmaI*reI_raw;
    gamma = prior_gamma_mean + sigmaGam*gamma_raw;
    delta = prior_delta_mean + sigmaDel*delta_raw;
    mui= bS + Z*gamma + Ht*delta + reI[subj_of_obs]; 
}
model {
    target += normal_lpdf(gamma_raw | 0, 1);
    target += normal_lpdf(delta_raw | 0, 1);
    target += normal_lpdf(reI_raw | 0, 1);
    target += normal_lpdf(sigmaI | prior_sigmaI_mean, prior_sigmaI_sd);
    target += normal_lpdf(sigmaGam | prior_sigmaGam_mean, prior_sigmaGam_sd);
    target += normal_lpdf(sigmaDel | prior_sigmaDel_mean, prior_sigmaDel_sd);
    target += binomial_logit_lpmf(y | nT, mui);
} 



