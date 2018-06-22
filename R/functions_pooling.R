

#' @title Measures of Variance Explained and Pooling/Shrinkage
#' @description Computes a form of R2 and Pooling Factor for a level of a hierarchical model
#' @param theta Matrix of the 'outcome' values for the level. Rows represent samples from the 
#' posterior distribution and columns represent units of the level.
#' @param epsilon Matrix of the 'error' values for the level, in the same format as \code{theta}.
#' @details This computes generalized versions of R^2 and pooling factor for levels of a hierarchical model, as described by Gelman and Pardoe (2006). 
#' @references Gelman, A and Pardoe, I. (2006). Bayesian Measures of Explained Variance and Pooling in Multilevel (Hierarchical) Models. \emph{Technometrics}, 48, 241-251.
#' @seealso \code{\link{exposure_pooling_factor}}
#' @export
#' 
compute_pooling_metrics <- function(theta, epsilon){
    c(prop_var = prop_var_explained(theta=theta, epsilon=epsilon),
      pooling_factor = pooling_factor(epsilon=epsilon))
}

#' @rdname compute_pooling_metrics
#' @export
prop_var_explained <- function(theta, epsilon){
    1 - mean(apply(epsilon, 1, var))/mean(apply(theta, 1, var))
}

#' @rdname compute_pooling_metrics
#' @export
pooling_factor <- function(epsilon){
    1 - var(apply(epsilon, 2, mean))/mean(apply(epsilon, 1, var))
}


#' @title Compute pooling factors for exposure model levels
#' @description Wrapper function for computing variance explained and pooling factor for a  fit exposure model.
#' @param stanfit The fitted model object from \code{\link{sample_exposure_model}}.
#' @param standata List containing model structure. Typically an \code{standata_outcome} object from \code{\link{create_standata_outcome}}.
#' @param level String indicating which level of the model to evaluate. Possible values are "cluster", "househould", "time", and "observation".
#' @details This is a wrapper around \code{\link{compute_pooling_metrics}} that extracts the necessary posterior samples from \code{stanfit}.
#' @seealso \code{\link{compute_pooling_metrics}} \code{\link{sample_exposure_model}} 
#' @export
exposure_pooling_factor <- function(stanfit, standata, level=c("observation", "household", "cluster", "group", "time")){
    if (level=="observation"){
        # Check for muWi......
        pm <- extract(stanfit, pars=c("etaK", "reH", "thetaG"))
        theta <- matrix(standata$w, nrow=nrow(pm$etaK), ncol=length(standata$w), byrow = TRUE)
        epsilon <- theta - (pm$etaK[, standata$cluster_of_obs] + pm$reH[, standata$hh_of_obs] + pm$thetaG %*% t(standata$Ht))
    } else if (level=="household"){
        # Check for muWi......
        pm <- extract(stanfit, pars=c("reH"))
        theta <- pm$reH
        epsilon <- pm$reH
    }  else if (level=="cluster"){
        # Check for muWi......
        pm <- extract(stanfit, pars=c("etaK","etaG"))
        theta <- pm$etaK
        epsilon <-  pm$etaK -  pm$etaG[, standata$group_of_cluster]
    } else if (level=="group"){
        # Check for muWi......
        pm <- extract(stanfit, pars=c("etaG"))
        theta <- pm$etaG - standata$prior_etaG_mean
        epsilon <-  pm$etaG
    } else if (level=="time"){
        # Check for muWi......
        pm <- extract(stanfit, pars=c("thetaG"))
        theta <- pm$thetaG - standata$prior_thetaG_mean
        epsilon <-  pm$thetaG
    } else {
        stop("not supported yet.")
    }
    compute_pooling_metrics(theta=theta, epsilon=epsilon)
}