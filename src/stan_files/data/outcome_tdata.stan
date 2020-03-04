//transformed data for outcome model
int<lower=1, upper=S> study_of_obs[N]; // which study is this subj in?
study_of_obs=study_of_unit[unit_of_obs];
