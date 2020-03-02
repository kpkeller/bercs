// Define transformed parameters (other than beta)
    vector[N] mui; // Mean of observations
    vector[n] reI; // Random effects by unitect
    vector[p] gamma; // Coefficients for covariates Z
    vector[timedf] delta; // Coefficients for time spline Ht
