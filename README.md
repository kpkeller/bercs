# bercs
Bayesian Exposure-Response Curves via STAN

This R package implements a two flexible hierarchical models. The *outcome model* is designed for estimating exposure-response curves, with particular focus on settings where data from multiple studies are being pooled. The *exposure model* is designed for modeling highly-variable, clustered, sparse longitudinal exposure measurements.

This package uses the Stan language and depends on the `rstan` package. To install `bercs`, use the following commands:  
```
devtools::install_github("jpkeller/bercs")
```


## Community guidelines

If you have a bug to report, are having technical issues, or want to recommend features, please open a [Github Issue](https://github.com/jpkeller/bercs/issues). 
