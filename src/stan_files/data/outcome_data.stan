// Data statements for outcome model
matrix[p==0 ? 0 : N,p] Z; // Matrix of covariates (not counting exposure, time, or intercept)
matrix[N,xdf] Mx; // Matrix of exposure values, transformed to basis
matrix[timedf==0 ? 0 : N,timedf] Mt;  // matrix of time spline values
int<lower=0> y[N]; // Case indicator
int<lower=0> nT[N]; // Time at risk
int<lower=1, upper=n> unit_of_obs[N]; // which subj is this obs?
int<lower=1, upper=S> study_of_subj[n]; // which study is this subj in?
