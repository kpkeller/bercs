

#' @title Measures of Variance Explained and Pooling/Shrinkage
#' @description Computes a form of R2 and Pooling Factor for a level of a hierarchical model
#' @param theta Matrix of the 'outcome' values for the level. Rows represent samples from the
#' posterior distribution and columns represent units of the level.
#' @param epsilon Matrix of the 'error' values for the level, in the same format as \code{theta}.
#' @details This computes generalized versions of R^2 and pooling factor for levels of a hierarchical model, as described by Gelman and Pardoe (2006).
#' @references Gelman, A and Pardoe, I. (2006). Bayesian Measures of Explained Variance and Pooling in Multilevel (Hierarchical) Models. \emph{Technometrics}, 48, 241-251.
#' @seealso \code{\link{exposure_pooling_factor}}
#' @export
compute_pooling_metrics <- function(theta, epsilon){
    c(prop_var = prop_var_explained(theta=theta, epsilon=epsilon),
      pooling_factor = pooling_factor(epsilon=epsilon))
}

#' @rdname compute_pooling_metrics
#' @export
#' @importFrom stats var
prop_var_explained <- function(theta, epsilon){
    1 - mean(apply(epsilon, 1, var))/mean(apply(theta, 1, var))
}

#' @rdname compute_pooling_metrics
#' @export
#' @importFrom stats var
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
#' @examples
#' data(casedataA)
#' data(casedataB)
#' outcome_combo_data <- create_standata_outcome(datalist=list(casedataA, casedataB),
#'                                               xdf=4,
#'                                               xfnargs=list(Boundary.knots=c(5, 200)))
#' outcome_combo_data <- add_priors(outcome_combo_data, sigmaI=c(0, 0.1))
#'
#' outcome_combo_mod_fit <- sample_outcome_model(outcome_combo_data,
#' B=2000,
#' cores=4)
#'
#' exposure_pooling_factor(stanfit = outcome_combo_mod_fit,
#' standata = outcome_combo_data,
#' level='observation')
#' @export
exposure_pooling_factor <- function(stanfit, standata, level=c("observation", "unit", "cluster", "group", "time")){
    if (level=="observation"){
        # Check for muWi......
        pm <- extract(stanfit, pars=c("reK", "reI", "theta", "etaG"))
        theta <- matrix(standata$w, nrow=nrow(pm$reI), ncol=length(standata$w), byrow = TRUE)
        epsilon <- theta - (pm$etaG[, standata$group_of_obs] + pm$reI[, standata$unit_of_obs] + pm$theta %*% t(standata$Ht))
        if ("reK" %in% names(pm)) epsilon <- epsilon -  pm$reK[, standata$cluster_of_obs]
    } else if (level=="unit"){
        pm <- extract(stanfit, pars=c("reI", "reK"))
        theta <- pm$reI
        if ("reK" %in% names(pm)) theta <- theta + pm$reK[, get_cluster_of_unit(standata)]
        epsilon <- pm$reI
    }  else if (level=="cluster"){
        pm <- extract(stanfit, pars=c("reK"))
        theta <- pm$reK
        epsilon <-  pm$reK
    } else if (level=="group"){
        pm <- extract(stanfit, pars=c("etaG"))
        theta <- pm$etaG
        epsilon <-  pm$etaG - standata$prior_etaG_mean
    } else if (level=="time"){
        pm <- extract(stanfit, pars=c("theta"))
        theta <- pm$theta
        epsilon <-  pm$theta - standata$prior_theta_mean
    } else {
        stop("not supported yet.")
    }
    compute_pooling_metrics(theta=theta, epsilon=epsilon)
}


get_cluster_of_unit <- function(standata){
    unit_of_obs <- standata$unit_of_obs
    inds <- !duplicated(unit_of_obs)
    cluster_of_unit <- standata$cluster_of_obs[inds]
    unit_of_obs <- unit_of_obs[inds]
    cluster_of_unit <- cluster_of_unit[order(unit_of_obs)]
    cluster_of_unit
}
