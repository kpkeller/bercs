/*
Exposure model

Note: No checks done for label consistency--user is responsible for this.
*/
data {
int<lower=1> G; // total number of treatment arms
int<lower=0> K; // total number of cluster
int<lower=0> n; // total number of units
int<lower=0> N; //total number of observations
int<lower=0> timedf; //degrees of freedom for time spline
// Labels
int<lower=1, upper=G> group_of_obs[N]; // which group is this obs in?
int<lower=0, upper=K> cluster_of_obs[N]; // which cluster is this obs in?
int<lower=1, upper=n> unit_of_obs[N]; // which unit is this obs for?
// Data
real w[N]; // observed values
matrix[timedf==0 ? 0 : N,timedf] Mt;  // matrix of time spline values
// Prior Specification
// Pass-in as data, to keep from needing to re-compile for each change
real prior_etaG_mean;
real prior_sigmaG_mean;
real<lower=0> prior_sigmaG_sd;
real prior_reK_mean;
real prior_sigmaK_mean;
real<lower=0> prior_sigmaK_sd;
real prior_reI_mean;
real prior_sigmaI_mean;
real<lower=0> prior_sigmaI_sd;
real prior_sigmaW_mean;
real<lower=0> prior_sigmaW_sd;
real prior_theta_mean;
real prior_sigmaTheta_mean;
real<lower=0> prior_sigmaTheta_sd;
}
parameters {
vector[G] etaG_raw; // transformed to etaG
vector[K] reK_raw; // transformed to etaK
vector[n] reI_raw; // transform to reI
vector[timedf] theta_raw;
real<lower=0> sigmaG; // SD of group means
real<lower=0> sigmaK[K==0 ? 0 : 1]; // SD of cluster RE
real<lower=0> sigmaI; // SD of unit RE
real<lower=0> sigmaW; // SD of measurements (instrument error)
real<lower=0> sigmaTheta[timedf==0 ? 0 : 1]; // SD of time trend coefs (instrument error)
}
transformed parameters {
vector[N] muW; // Mean of observations
vector[G] etaG;
vector[K] reK;
vector[n] reI;
vector[timedf] theta;
etaG = prior_etaG_mean + sigmaG * etaG_raw;
reI = prior_reI_mean + sigmaI * reI_raw;
muW = etaG[group_of_obs] + reI[unit_of_obs];
if (K > 0 ){
    reK = prior_reK_mean + sigmaK[1] * reK_raw;
    muW += reK[cluster_of_obs];
}
if (timedf > 0){
    theta = prior_theta_mean + sigmaTheta[1]*theta_raw;
    muW += Mt*theta;
}
}
model {
if (timedf > 0){
    target += normal_lpdf(theta_raw | 0, 1);
    target += normal_lpdf(sigmaTheta | prior_sigmaTheta_mean, prior_sigmaTheta_sd);
}
if (K > 0 ){
target += normal_lpdf(reK_raw | 0, 1);
target += normal_lpdf(sigmaK | prior_sigmaK_mean, prior_sigmaK_sd);
}
target += normal_lpdf(sigmaI | prior_sigmaI_mean, prior_sigmaI_sd);
target += normal_lpdf(reI_raw | 0, 1);
target += normal_lpdf(sigmaG | prior_sigmaG_mean, prior_sigmaG_sd);
target += normal_lpdf(etaG_raw | 0, 1);
target += normal_lpdf(sigmaW | prior_sigmaW_mean, prior_sigmaW_sd);
target += normal_lpdf(w | muW, sigmaW);
}
generated quantities {
real<lower=0> sigma2I;
real<lower=0> sigma2W;
real<lower=0> sigma2G;
real<lower=0> sigma2K[K==0 ? 0 : 1]; // Var of cluster RE
sigma2I=square(sigmaI);
sigma2W=square(sigmaW);
sigma2G=square(sigmaG);
if (K > 0){
    sigma2K=square(sigmaK);
}
}
