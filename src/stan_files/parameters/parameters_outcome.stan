// Outcome model parameters
// other than beta

real<lower=0> sigmaI;
real<lower=0> sigmaGamma[p==0 ? 0 : 1] ;
real<lower=0> sigmaDelta[timedf==0 ? 0 : 1] ;
vector[S] bS; // Study-level intercept
vector[p] gamma_raw; // Coefficients for covariates Z
vector[timedf] delta_raw; // Coefficients for time spline Ht
vector[n] reI_raw;
