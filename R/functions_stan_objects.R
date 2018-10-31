



# Wrapper function for adding prior values to standata object fed in to STAN.
##' @title Add Hyperparameter Values to Standata Object
##' @description Sets prior distribution hyperparameters using defaults or provided values
##' @param standata Either a \code{standata_exposure} or \code{standata_outcome} object. Based on its class, \code{add_priors} calls \code{add_priors_exposure_model} or \code{add_priors_exposure_model}.
##' @param ... Arguments passed to \code{add_priors_exposure_model} or \code{add_priors_outcome_model}. Unlisted arguments are ignorred.
##' @seealso \code{\link{create_standata_exposure}}
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


# Function for adding prior values to standataect for Exposure model
# Numeric pairs are (mean, sd) for variables
##' @rdname add_priors
#' @param etaG Prior mean and standard deviation for the group means parameter `etaG`.
#' @param sigmaG Prior mean and standard deviation for the group standard deviation parameter `sigmaG`.
#' @param reK Prior mean for cluster-level random effects.
#' @param reH Prior mean for household-level random effects.
#' @param sigmaK Prior mean and standard deviation for the standard deviation of the cluster random effects, `sigmaK`.
#' @param sigmaH Prior mean and standard deviation for the standard deviation of the household random effects, `sigmaH`.
#' @param sigmaW Prior mean and standard deviation for the residual standard deviation parameter `sigmaW`.
#' @param theta Prior mean and standard deviation for the time coefficients `thetaG`.
#' @param sigmaTheta Prior mean and standard deviation for the standard deviation parameter for time trends, `sigmaTheta`.
##' @details For \code{add_stan_priors_exposure_model}, the prior distributions are normal distributions parameterixed by (mean, sd).
##' @export
add_priors_exposure_model <- function(standata, etaG=0, sigmaG=c(0, 1), reK=0, sigmaK=c(0, 1), reH=0,sigmaH=c(0, 1), sigmaW=c(0, 1),  theta=0, sigmaTheta=c(0, 1), ...) {

    standata$prior_etaG_mean <- etaG[1]
    standata$prior_sigmaG_mean <- sigmaG[[1]]
    standata$prior_sigmaG_sd <- sigmaG[[2]]

    standata$prior_sigmaW_mean <- sigmaW[[1]]
    standata$prior_sigmaW_sd <- sigmaW[[2]]

    standata$prior_reH_mean <- reK[1]
    standata$prior_sigmaH_mean <- sigmaH[[1]]
    standata$prior_sigmaH_sd <- sigmaH[[2]]

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
##' @param sigmaGam Prior mean and standard deviation for the standard deviation of covariate coefficient parameter.
##' @param delta Prior mean fro the time coefficient(s).
##' @param sigmaDel Prior mean and standard deviation for the standard deviation of time coefficient parameter.
##' @param sigmaI Prior mean and standard deviation for the standard deviation of the subject-level random effect.
##' @param sigmaK Prior mean and standard deviation for the standard deviation of the cluster-level random effect.
##' @param beta_nu Prior value for the LKJ prior on the correlation between exposure coefficients from different studies.
##' @export
add_priors_outcome_model <- function(standata, beta=0, sigmaBeta=c(0, 1), gamma=0, sigmaGam=c(0 ,1), delta=0, sigmaDel=c(0, 1), sigmaI=c(0, 1), sigmaK=c(0, 1), beta_nu=1, ...) {

    standata$prior_beta_mean <- beta
    standata$prior_sigmaBeta_mean <- sigmaBeta[[1]]
    standata$prior_sigmaBeta_sd <- sigmaBeta[[2]]
    standata$beta_nu <- beta_nu

    standata$prior_gamma_mean <- gamma
    standata$prior_sigmaGam_mean <- sigmaGam[[1]]
    standata$prior_sigmaGam_sd <- sigmaGam[[2]]

    standata$prior_delta_mean <- delta
    standata$prior_sigmaDel_mean <- sigmaDel[[1]]
    standata$prior_sigmaDel_sd <- sigmaDel[[2]]

    standata$prior_sigmaI_mean <- sigmaI[[1]]
    standata$prior_sigmaI_sd <- sigmaI[[2]]

    standata$prior_sigmaK_mean <- sigmaK[[1]]
    standata$prior_sigmaK_sd <- sigmaK[[2]]

    standata
}



