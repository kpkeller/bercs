# bercs
Bayesian Exposure-Response Curves via STAN

This R package implements a flexible hierarchical model for estimating exposure-response curves.

This package uses the Stan language and depends on the `rstan` package. To install `bercs`, use the following commands:  
```
Sys.setenv(USE_CXX14 = 1)
devtools::install_github("jpkeller/bercs")
```