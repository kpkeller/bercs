// Outcome model dimensions

int<lower=1> S; //total number of studies
int<lower=1> N; //total number of observations
int<lower=1> n; //total number of units
int<lower=0> xdf; // number of basis functions for exposure (e.g. PM)
int<lower=0> timedf; // degrees of freedom in time spline
int<lower=0> p; // number of covariates -- do NOT include intercept here
