
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bercs

Bayesian Exposure-Response Curves via STAN

This R package implements a two flexible hierarchical models. The
**exposure model** is designed for modeling highly-variable, clustered,
sparse longitudinal exposure measurements. The **outcome model** is
designed for estimating exposure-response curves, with particular focus
on settings where data from multiple studies are being pooled.

This package uses the Stan language and depends on the `rstan` package.
To install `bercs`, use the following commands:

    devtools::install_github("jpkeller/bercs")

## Exposure Model

The exposure model is a hierarchical model for concentrations as a
function of group means, cluster random effects, unit (e.g. subject or
household) random effects, and temporal spline. The general procedure
for fitting the model is the following:

1.  Create the data object using `create_standata_exposure()`. This
    function creates a list with the named components used in model
    fitting.

<!-- end list -->

``` r
# Create simulated data with:
#  3 groups
#  each with 6 units
#  that each have 4 observations
#  at times 1, 2, 3, 4
exposure_standata <- create_standata_exposure(group=rep(1:3, each=24),
                                              conc=rnorm(n=72,mean=rep(c(0, 2, 4), each=23)),
                                              unit_id=rep(1:18, each=4),
                                              time=rep(1:4, times=18))
```

2.  If using a time spline, add the corresponding matrix of spline
    values to the model.

<!-- end list -->

``` r
# Add natural splines with 3 df
Mt <- create_spline_matrix(x=exposure_standata$time,
                    df=3,
                    fn="ns")
exposure_standata <- add_spline_time(exposure_standata, Mt=Mt)
```

3.  Add the hyperparameters for the prior distributions.

<!-- end list -->

``` r
# Use defaults for all except sigmaI
exposure_standata <- add_priors(exposure_standata,
                                sigmaI=c(0, 0.5))
```

4.  Sample from the posterior distribution of parameters using STAN.

<!-- end list -->

``` r
exposure_mod_fit <- sample_exposure_model(exposure_standata,
                                          B=1000,
                                          chains=4)
```

``` r
print(exposure_mod_fit,
      pars=c("muW", "reI_raw", "reI", "etaG_raw", "theta_raw"),
      include=FALSE)
```

    ## Inference for Stan model: exposure_model.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##                  mean se_mean   sd    2.5%     25%     50%     75%   97.5%
    ## sigmaG           1.84    0.01 0.46    1.12    1.50    1.77    2.12    2.91
    ## sigmaI           0.63    0.01 0.24    0.12    0.49    0.64    0.79    1.09
    ## sigmaW           1.20    0.00 0.13    0.99    1.11    1.19    1.28    1.48
    ## sigmaTheta[1]    0.35    0.01 0.33    0.01    0.12    0.25    0.48    1.24
    ## etaG[1]          0.03    0.01 0.37   -0.70   -0.21    0.03    0.27    0.79
    ## etaG[2]          1.95    0.01 0.38    1.17    1.71    1.96    2.21    2.67
    ## etaG[3]          3.38    0.01 0.38    2.60    3.14    3.39    3.64    4.10
    ## theta[1]         0.04    0.00 0.29   -0.55   -0.09    0.01    0.15    0.70
    ## theta[2]        -0.05    0.00 0.32   -0.80   -0.17   -0.01    0.09    0.63
    ## theta[3]         0.05    0.00 0.23   -0.39   -0.07    0.02    0.15    0.57
    ## sigma2I          0.46    0.01 0.31    0.01    0.24    0.40    0.62    1.19
    ## sigma2W          1.46    0.01 0.31    0.98    1.24    1.42    1.65    2.18
    ## sigma2G          3.60    0.04 1.89    1.25    2.25    3.13    4.49    8.48
    ## lp__          -158.37    0.18 5.37 -170.17 -161.62 -157.95 -154.48 -149.23
    ##               n_eff Rhat
    ## sigmaG         2275    1
    ## sigmaI         1290    1
    ## sigmaW         1948    1
    ## sigmaTheta[1]  3424    1
    ## etaG[1]        2536    1
    ## etaG[2]        3689    1
    ## etaG[3]        3654    1
    ## theta[1]       4836    1
    ## theta[2]       5016    1
    ## theta[3]       4996    1
    ## sigma2I        1730    1
    ## sigma2W        1910    1
    ## sigma2G        2367    1
    ## lp__            858    1
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Jun 18 16:14:11 2020.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

5.  Compute long-term means and plot them.

<!-- end list -->

``` r
fitted_means <- compute_fitted_mean(stanfit=exposure_mod_fit,
                                    standata=exposure_standata)
summary(fitted_means)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## -0.5347  0.4131  1.8332  1.8291  3.5621  3.9226

``` r
# Plots the fitted values (colored points) on top of original data (black outlines).
plot_exposure_means_bytime(stanfit=exposure_mod_fit,
                           standata=exposure_standata)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Exposure-Response (Outcome) Model

The outcome model is a hierarchical model for estimating an
exposure-response function. It can accommodate data from multiple
studies, each with their own measured covariates, temporal trend, and
overall risk mean. Currently, the model is designed only for binary
outcomes (i.e. logistic regression). In addition to the study-specific
time trends and covariate effects, subject-level random effects can be
included. The general procedure for fitting the model is the following:

1.  Create the data object using `create_standata_outcome()`. This
    function creates a list with the named components used in model
    fitting.

<!-- end list -->

``` r
# Dataset A has a linear exposure-response relationship
# across the exposure range of 5-50
data(casedataA)
# Dataset B has no exposure-response relationship
# and an exposure range of 50-200
data(casedataB)
outcome_combo_data <- create_standata_outcome(datalist=list(casedataA, casedataB),
                                              xdf=4,
                                              xfnargs=list(Boundary.knots=c(5, 200),
                                                           knots=c(25)))
```

2.  Add the hyperparameters for the prior distributions.

<!-- end list -->

``` r
# Use defaults for all except sigmaI
outcome_combo_data <- add_priors(outcome_combo_data,
                                 sigmaI=c(0, 0.1))
```

3.  Sample from the posterior distribution of parameters using STAN.

<!-- end list -->

``` r
outcome_combo_mod_fit <- sample_outcome_model(outcome_combo_data,
                                         B=2000)
```

``` r
print(outcome_combo_mod_fit, pars=c("reI_raw", "reI","mui", "beta_raw"), include=FALSE)
```

    ## Inference for Stan model: outcome_model.
    ## 4 chains, each with iter=4000; warmup=2000; thin=1; 
    ## post-warmup draws per chain=2000, total post-warmup draws=8000.
    ## 
    ##              mean se_mean    sd     2.5%     25%     50%     75%   97.5% n_eff
    ## sigmaI       0.08    0.00  0.06     0.00    0.03    0.07    0.11    0.22  6789
    ## bS[1]       -1.65    0.01  0.38    -2.45   -1.90   -1.63   -1.38   -0.98  3673
    ## bS[2]       -2.93    0.01  0.76    -4.45   -3.45   -2.94   -2.39   -1.52  3819
    ## sigmaBeta    1.21    0.01  0.59     0.15    0.79    1.17    1.59    2.48  3626
    ## beta[1]      1.91    0.02  1.11     0.01    1.10    1.85    2.65    4.27  3449
    ## beta[2]     -0.43    0.01  1.01    -2.78   -0.97   -0.26    0.20    1.29  6622
    ## beta[3]     -0.23    0.01  0.92    -2.20   -0.75   -0.17    0.27    1.60  7352
    ## beta[4]      0.38    0.01  0.85    -1.17   -0.14    0.28    0.86    2.27  7659
    ## lp__      -971.91    0.30 16.06 -1003.99 -982.44 -971.85 -961.07 -940.77  2861
    ##           Rhat
    ## sigmaI       1
    ## bS[1]        1
    ## bS[2]        1
    ## sigmaBeta    1
    ## beta[1]      1
    ## beta[2]      1
    ## beta[3]      1
    ## beta[4]      1
    ## lp__         1
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Jun 18 16:14:47 2020.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

5.  Plot the exposure-response curve

<!-- end list -->

``` r
fitted_ERC <- get_fitted_ERC(standata=outcome_combo_data,
                             stanfit=outcome_combo_mod_fit,
                             exprange=c(5,200))
plot_fitted_ERC(fitted_ERC) + scale_y_log10()
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## Community guidelines

If you have a bug to report, are having technical issues, or want to
recommend features, please open a [Github
Issue](https://github.com/jpkeller/bercs/issues).
