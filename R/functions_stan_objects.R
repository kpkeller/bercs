



# Wrapper function for adding prior values to standata object fed in to STAN.
##' @title Add Hyperparameter Values to Standata Object
##' @description Sets prior distribution hyperparameters using defaults or provided values
##' @param standata Either a \code{standata_exposure} or \code{standata_outcome} object.
##' @param ... Arguments passed to \code{add_priors_exposure_model} or \code{add_priors_outcome_model}. Unlisted arguments are ignored.
##' @details Based on its class, \code{add_priors} calls \code{add_priors_exposure_model} or \code{add_priors_exposure_model}. For the variance parameters, the prior distributions are half-normal distributions parameterized by (mean, sd).
##' @seealso \code{\link{create_standata_exposure}}, \code{\link{create_standata_outcome}}, \code{\link{sample_exposure_model}}, \code{\link{sample_outcome_model}}
##' @export
add_priors  <- function(standata, ...){
    if (inherits(standata, "standata_exposure")){
        out <- add_priors_exposure_model(standata, ...)
    } else if (inherits(standata, "standata_outcome")){
        out <- add_priors_outcome_model(standata, ...)
    } else {
        stop("`standata` is not of a supported class.")
    }
    out
}

# Function for adding prior values to standata for Exposure model
# Numeric pairs are (mean, sd) for variables
##' @rdname add_priors
#' @param etaG Prior mean and standard deviation for the group means parameter `etaG`.
#' @param sigmaG Prior mean and standard deviation for the group standard deviation parameter `sigmaG`.
#' @param reK Prior mean for cluster-level random effects.
#' @param reI Prior mean for household-level random effects.
#' @param sigmaK Prior mean and standard deviation for the standard deviation of the cluster random effects, `sigmaK`.
#' @param sigmaI Prior mean and standard deviation for the standard deviation of the household random effects, `sigmaI`.
#' @param sigmaW Prior mean and standard deviation for the residual standard deviation parameter `sigmaW`.
#' @param theta Prior mean and standard deviation for the time coefficients `thetaG`.
#' @param sigmaTheta Prior mean and standard deviation for the standard deviation parameter for time trends, `sigmaTheta`.
##' @export
##' @examples
##' # Create simulated data
##' exp_data <- create_standata_exposure(group=rep(1, 10),
##'                                      conc=rnorm(10),
##'                                      unit_id=rep(0:1, 5),
##'                                      time=runif(10))
##' # Add comnbination of default and custom prior
##' exp_data <- add_priors(exp_data,
##'                        sigmaI=c(0, 0.1),
##'                        sigmaK=c(0, 2))
add_priors_exposure_model <- function(standata, etaG=0, sigmaG=c(0, 1), reK=0, sigmaK=c(0, 1), reI=0,sigmaI=c(0, 1), sigmaW=c(0, 1),  theta=0, sigmaTheta=c(0, 1), ...) {

    standata$prior_etaG_mean <- etaG[1]
    standata$prior_sigmaG_mean <- sigmaG[[1]]
    standata$prior_sigmaG_sd <- sigmaG[[2]]

    standata$prior_sigmaW_mean <- sigmaW[[1]]
    standata$prior_sigmaW_sd <- sigmaW[[2]]

    standata$prior_reI_mean <- reK[1]
    standata$prior_sigmaI_mean <- sigmaI[[1]]
    standata$prior_sigmaI_sd <- sigmaI[[2]]

    standata$prior_reK_mean <- reK[1]
    standata$prior_sigmaK_mean <- sigmaK[[1]]
    standata$prior_sigmaK_sd <- sigmaK[[2]]

    standata$prior_theta_mean <- theta[1]
    standata$prior_sigmaTheta_mean <- sigmaTheta[[1]]
    standata$prior_sigmaTheta_sd <- sigmaTheta[[2]]

    standata
}


##' @rdname add_priors
##' @param beta Prior mean for the exposure coefficient(s).
##' @param sigmaBeta Prior mean and standard deviation for the standard deviation of exposure coefficient parameter.
##' @param gamma Prior mean for the covariate coefficient(s).
##' @param sigmaGamma Prior mean and standard deviation for the standard deviation of covariate coefficient parameter.
##' @param delta Prior mean for the time coefficient(s).
##' @param sigmaDelta Prior mean and standard deviation for the standard deviation of time coefficient parameter.
##' @param sigmaI Prior mean and standard deviation for the standard deviation of the subject-level random effect.
##' @param beta_nu Prior value for the LKJ prior on the correlation between exposure coefficients from different studies.
##' @param sigma_y Prior mean and standard deviation for continuous outcome y.
##' @export
add_priors_outcome_model <- function(standata, beta=0, sigmaBeta=c(0, 1), gamma=0, sigmaGamma=c(0 ,1), delta=0, sigmaDelta=c(0, 1), sigmaI=c(0, 1), beta_nu=1, sigma_y=c(0, 1), ...) {

    standata$prior_beta_mean <- beta
    standata$prior_sigmaBeta_mean <- sigmaBeta[[1]]
    standata$prior_sigmaBeta_sd <- sigmaBeta[[2]]
    standata$beta_nu <- beta_nu

    standata$prior_gamma_mean <- gamma
    standata$prior_sigmaGamma_mean <- sigmaGamma[[1]]
    standata$prior_sigmaGamma_sd <- sigmaGamma[[2]]

    standata$prior_delta_mean <- delta
    standata$prior_sigmaDelta_mean <- sigmaDelta[[1]]
    standata$prior_sigmaDelta_sd <- sigmaDelta[[2]]

    standata$prior_sigmaI_mean <- sigmaI[[1]]
    standata$prior_sigmaI_sd <- sigmaI[[2]]

    standata$prior_sigma_y_mean <- sigma_y[[1]]
    standata$prior_sigma_y_sd <- sigma_y[[2]]

    standata
}



