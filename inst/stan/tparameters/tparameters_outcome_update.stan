// Update to transformed parameters
// other than beta's
reI = sigmaI*reI_raw;
mui= bS[study_of_obs] + reI[unit_of_obs];
if (p >0 ){
  gamma = prior_gamma_mean + sigmaGamma[1]*gamma_raw;
  mui += Z*gamma;
}
if (timedf > 0){
  delta = prior_delta_mean + sigmaDelta[1]*delta_raw;
  mui += Mt*delta;
}
