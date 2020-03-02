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
target += binomial_logit_lpmf(y | nT, mui);
