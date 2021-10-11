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

