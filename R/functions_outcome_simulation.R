# Outcome model


##' @title Create Skeleton for Outcome Simulations
##' @description These functions construct lists containing the structure and data for simulating exposure and outcome data.
##' @param design String indicating trial design. Currently only 'parallel' is supported.
##' @param ... Additional arguments passed to the design-specific functions.
##' @details These functions create the blank structure of groups, clusters, households, and subject.
##' Once parameters are set via \code{\link{expsim_update_parameter}} and related functions, then data can be sampled via 
##' \code{\link{expsim_sample_observations}}
##' and posterior parameter estimates obtained from \code{\link{sample_exposure_model}}.
##' 
##' For creating a \code{standata_outcome} object from an existing 'wide' format data frame, see \code{\link{create_standata_outcome}}.
##' @return A list containing two sublists: \code{structure}, which contains settings and parameters for generating data, and \code{standata} which contains the study data in a format for sampling via STAN.
#' @family outcome simulation functions
##' @seealso \code{\link{sample_outcome_model}}, \code{\link{create_exposure_simulation_skeleton}}
##' @export
create_outcome_simulation_skeleton <- function(design="parallel",...){
    if (design=="parallel"){
        obj <- create_outcome_simulation_skeleton_parallel(...)
    } else if (design=="crossover"){
        stop("Crossover design not yet implemented.")
        # obj <- expsim_create_study_skeleton_crossover(...)
    } else {
        stop(paste0("Design = ", design, " is not supported. Please use 'parallel'."))
    }
    obj
}



# Need to have clear method for converting x to X
# and error checking for x
##' @rdname create_outcome_simulation_skeleton
##' @param nstudies Number of studies.
##' @param nclusters Number of clusters per study. If not of length \code{nstudies}, the same value will be used for all studies.
##' @param nsubjects Number of subjects.
##' @param nobs Number of observations per subject. If length 1, value is repeated for all subjects. If length \code{sum(nclusters)}, then each element is repeated for all subjects within the corresponding cluster. 
##' @param study_of_cluster optional vector of positive integers that provide the study number of each cluster. Should have length \code{sum(nclusters)}.
##' @param cluster_of_subj optional vector that provides the cluster number of each subject. Should have length \code{sum(nsubjects)}.
##' @param subj_of_obs optional vector that provides the subject number of each observation. Should have length \code{sum(nobs)}.
##' @param beta0 Study-level intercept, on the logit scale. Defaults to 0 and can be updated later via \code{\link{outsim_update_covariate}}.
##' @param x Exposure concentration values. Defaults to 0 and can be updated later via \code{\link{outsim_update_covariate}}.
##' @param nT optional time-at-risk value. See \code{\link{outsim_update_atrisk}}
##' @param times optional times for observations. Can be set later via \code{\link{sim_update_times}}.
##' @param xfn Exposure-response function. Can be set later via \code{\link{outsim_update_parameter}}. 
##' @param timefn Function relating time to mean outcome. Can be set later via \code{\link{outsim_update_parameter}}.
##' @param verbose should messages be printed.
##' @export
create_outcome_simulation_skeleton_parallel <- function(nstudies=1,
                                                  nclusters=1,
                                                  nsubjects=1,
                                                  nobs=1,
                                                  study_of_cluster,
                                                  cluster_of_subj,
                                                  subj_of_obs,
                                                  beta0=0,
                                                  x=0,
                                                  nT=NA,
                                                  times=NA,
                                                  verbose=TRUE,
                                                  xfn=function(x) {x},
                                                  timefn=function(t) {rep(0, length(t))}){
    # Check input
    if(!is.null(nclusters) && !check_all_nonnegative_value(nclusters)) stop("'nclusters' should be a non-negative integer and must have at least one positive element.")
    if(!check_all_nonnegative_value(nsubjects)) stop("'nsubjects' should be a non-negative integer and must have at least one positive element.")
    if(!check_all_postive_value(nstudies)) stop("'nstudies' should be a positive integer.")
    
    # Process nclusters
    if(!length(nclusters) %in% c(1, nstudies)) stop("'nclusters' must have length 1 or 'nstudies'.")
    if (length(nclusters)==1) {
        nclusters <- rep(nclusters, times=nstudies)
        if(verbose) {
            message("Only one value of 'nclusters' provided. Applying value to all studies")
        }
    }
    K <- sum(nclusters)
    
    # Process nsubjects
    if(!length(nsubjects) %in% c(1, K)) stop("'nsubjects' must have length 1 or 'sum(nclusters)'.")
    if (length(nsubjects)==1) {
        nsubjects <- rep(nsubjects, times=K)
        if(verbose) {
            message("Only one value of 'nsubjects' provided. Applying value to all clusters.")
        }
    }
    n <- sum(nsubjects)
    
    # Process nobs
    if(!length(nobs) %in% c(1, K, n)) stop("'nobs' must have length 1, 'sum(nclusters)', or 'sum(nsubjects)'.")
    if (length(nobs)==1) {
        nobs <- rep(nobs, times=n)
        if(verbose) {
            message("Only one value of 'nobs' provided. Applying value to all subjects.")
        }
    } else if (length(nobs)==K) {
        nobs <- rep(nobs, times=nsubjects)
        if(verbose) {
            message("One value of 'nobs' provided for each cluster. Applying each value to all subjects within corresponding cluster.")
        }
    }
    N <- sum(nobs)
    
    # Defaults
    if (missing(study_of_cluster)) study_of_cluster <- rep(1:nstudies, times=nclusters)
    if (missing(cluster_of_subj)) cluster_of_subj <- rep(1:K, times=nsubjects)
    if (missing(subj_of_obs)) subj_of_obs <- rep(1:n, times=nobs)
    cluster_of_obs <- cluster_of_subj[subj_of_obs]
    study_of_subj <- study_of_cluster[cluster_of_subj]
    study_of_obs <- study_of_cluster[cluster_of_obs]
    
    # Defaults
    if (is.na(nT)) nT <- rep(1, N)
    if (length(x)==1) x <- rep(x, N)
    if (length(times)==1) times <- rep(times, N)
    if (length(beta0)==1) beta0 <- rep(beta0, nstudies)
    
    study_structure <- list(nstudies=nstudies,
                            nclusters=nclusters,
                            nsubjects=nsubjects,
                            nobs=nobs,
                            study_of_cluster=study_of_cluster,
                            study_of_obs=study_of_obs,
                            cluster_of_subj=cluster_of_subj,
                            subj_of_obs=subj_of_obs,
                            cluster_of_obs=cluster_of_obs,
                            beta0=beta0,
                            gamma=NA,
                            sigI=NA,
                            reI=rep(NA, n),
                            times=times,
                            x=x,
                            xfn=xfn,
                            timefn=timefn,
                            logitmean=rep(NA, N))# , # mean for generating data, updated before return
    study_standata <- list(S=nstudies,
                           K=K,
                           n=n,
                           N=N,
                           Ht=matrix(0, 0, 0),
                           Ht_attributes=NA,
                           timedf=0,
                           y=rep(NA, N),
                           nT=nT,
                           x=x,
                           Hx=NA,
                           Hx_attributes=NA,
                           xdf=0,
                           Z=matrix(0, 0, 0),
                           p=0,
                           study_of_cluster=study_of_cluster,
                           cluster_of_subj=cluster_of_subj,
                           subj_of_obs=subj_of_obs,
                           study_of_subj=study_of_subj,
                           study_of_obs=study_of_obs,
                           cluster_of_obs=cluster_of_obs)
    
    class(study_standata) <- "standata_outcome"
    
    obj <- list(structure=study_structure,
                standata=study_standata)
    class(obj) <- "outsim"
    
    obj <- outsim_update_logitmean(obj)
    obj
}



##' @title Sample Simulated Outcome Observations
##' @description Draws a sample of outcome observations using provided model parameters
##' @param obj the outcome simulation object for which to sample
##' @family outcome simulation functions
##' @seealso \code{\link{expsim_sample_observations}}
##' @export
##' @importFrom stats rbinom
outsim_sample_observations <- function(obj){
    if (!inherits(obj, "outsim")) stop("'obj' must be of class 'outsim'.")
    obj <- outsim_update_logitmean(obj)
    obj$standata$y=stats::rbinom(n=obj$standata$N, size=obj$standata$nT, prob=expit(obj$structure$logitmean))
    obj
}

# Internal function to update the logitmean
outsim_update_logitmean <- function(obj){
    if (!inherits(obj, "outsim")) stop("'obj' must be of class 'outsim'.")
    obj$structure$logitmean <- with(obj$structure, beta0[study_of_obs] + xfn(x) + timefn(times))
    if(!any(is.na(obj$structure$gamma)) && !any(is.na(obj$standata$Z))){
        obj$structure$logitmean <- with(obj$structure, logitmean + obj$standata$Z %*% gamma)
    }
    if (!any(is.na(obj$structure$reI))){
        obj$structure$logitmean <- with(obj$structure, logitmean + reI[subj_of_obs])
    }
    obj
}



##' @title Update Outcome model parameters
##' @description Sets or overwrites outcome model parameters
#' @param obj 'outsim' object to update
#' @param param character string giving name of parameter to update. If not provided, parameter name is taken from \code{level} and \code{type}
#' @param level character string of model level being updated. See details.
#' @param type character string giving type of parameter being updated. See details.
#' @param val value(s) to which parameter should be set
#' @param draw if \code{TRUE}, then the parameter value is sampled instead of being set to \code{val}
#' @details Only specific combinations of \code{level} and \code{type} are allowed. They are:
#' \itemize{
#' \item \code{"group"}: \code{"mean"} and \code{"sd"} set 'etaG' and 'sigG', respectively
#' \item \code{"cluster"}: \code{"mean"} and \code{"sd"} set 'etaK' and 'sigK', respectively
#' }
#' @family outcome simulation functions
#' @seealso \code{\link{expsim_update_parameter}}
##' @export
outsim_update_parameter <- function(obj, param, level=c("study", "subject", "time", "covariate", "exposure"), type=c("mean", "sd", "re", "coef"), val=NULL, draw=is.null(val)){
    if (!inherits(obj, c("outsim"))) stop("'obj' must be of class 'outsim'.")

    if (!missing(param)){
        param_list <- c("beta0", "sigI", "reI", "timefn", "gamma", "xfn") #c("etaG", "etaK", "sigG", "sigK", "sigH", "reH",  "sigW", "corHW", "timefn")
        param_check <- param %in% param_list
        if (!param_check) stop('"param" must be one of  c("beta0", "sigI", "reI", "timefn","gamma", "xfn")')
    } else {
        level <- match.arg(level)
        type <- match.arg(type)
        param <- switch(paste0(level, "_", type),
                        study_mean="beta0",
                        subject_sd="sigI",
                        subject_re="reI",
                        time_mean="timefn",
                        covariate_coef="gamma",
                        exposure_mean="xfn",
                        NA)
        if (is.na(param)) stop(paste0(level, "_", type, " is not a supported parameter."))
    }
    
    # Need to do some error/length checking.
    current_length <- length(obj$structure[[param]])
    if (draw){
        if(!is.null(val)) warning("'val' provided but 'draw=TRUE'. Ignoring provided value of 'val'.")
        
        if (!param %in% c("reI")) {
            stop(paste0(param, " should be set, not drawn from a distribution."))
        }
        
        # Need to check that mean/sd has been set, first.
        cur_mean <- switch(param,
                           reI=0)
        cur_sd <- switch(param,
                         reI=obj$structure$sigI)
        cur_n <- switch(param,
                        reI=obj$standata$n)
        val <- rnorm(cur_n, cur_mean, cur_sd)
    }
    if (length(val)!=current_length) {
        if (length(val)==1){
            val <- rep(val, current_length)
        } else {
            stop(paste0("'val' should have length ", current_length, "."))
        }
    }
    obj$structure[[param]] <- val
    
    obj
}


##' @rdname outsim_update_parameter
##' @param nT time at risk value to set
##' @export
outsim_update_atrisk <- function(obj, nT=1){
    if (!inherits(obj, c("outsim"))) stop("'obj' must be of class 'outsim'.")
    if (length(nT)==1) nT <- rep(nT, obj$N)
    if (length(nT)!=obj$N) stop("Incorrect length for nT.")
    obj$nT <- nT    
    obj
}





##' @rdname outsim_update_parameter
##' @param covariate character string indicating which covariate to update. currently only 'x' and 'Z'.
##' @param ncol number of columns in covariate matrix. only used if \code{draw=TRUE}.
##' @param verbose Set to \code{FALSE} to suppress messages.
##' @param mean Mean of normal distribution from which to draw.
##' @param sd SD of normal dsitirbtion from which to draw.
##' @export
##' @importFrom stats rnorm
outsim_update_covariate <- function(obj, covariate=c("x","Z"), ncol=1, val=NULL, draw=is.null(val), verbose=TRUE, mean=0, sd=1){
    covariate <- match.arg(covariate)
    if (covariate=="x" && ncol!=1){
        stop("'ncol' must equal 1 when 'covariate' is 'x'.")
    }
    if (draw){
        if(!is.null(val)) warning("'val' is provided, but 'draw=TRUE'. Ignoring provided value of 'val'.")

        if (verbose){
            message(paste0("Drawing N(", mean, ",", sd, ") sample for covariate ", covariate, "."))
        }
        valdim <- c(obj$standata$N, ncol)
        val <- matrix(stats::rnorm(valdim[1]*valdim[2], mean=mean, sd=sd), nrow=valdim[1], ncol=valdim[2])
    } else {
        valdim <- ifelse(is.vector(val), c(length(val), 1), dim(val))
    }
    covariate_param <- switch(covariate,
                              Z="gamma",
                              NA)
    if (!is.null(obj$structure[[covariate_param]])){
        if (length(obj$structure[[covariate_param]])==1){
            obj$structure[[covariate_param]] <- rep(obj$structure[[covariate_param]], valdim[2])
        } else {
            obj$structure[[covariate_param]] <- obj$structure[[covariate_param]][1:valdim[2]]
            if (draw){
                warning(paste0("Truncating model parameter ",covariate_param, " to the first ncol=", ncol, " elements."))
            } else {
                warning(paste0("Truncating model parameter ",covariate_param, " to the first ncol(val)=", valdim[2], " elements."))
            }
        }
    }
    if (covariate=="Z") obj$standata$p <- valdim[2]
    
    obj$structure[[covariate]] <- val
    obj$standata[[covariate]] <- val
    obj
}

##' @rdname outsim_update_parameter
##' @export
outsim_update_set_parameter <- function(obj, param=c("beta", "gamma"), val=NULL){
    param <- match.arg(param)
    
    if (is.null(val)){
        stop("'val' must be provided.")
    }
    
    # Matrix corresponding to this parameer
    param_mx <- switch(param, 
                       beta="X",
                       gamma="Z",
                       beta0=NA)
    
    # Length of val should match dimensions of corresponding matrix
    # Repeat val if needed
    if (is.na(param_mx)){
        # if param=beta0, then val should have length 1.
        val_length <- 1
    } else if (!is.null(obj$standata[[param_mx]])) {
        # if matrix exists, val should have length 1, or number of columns of matrix
        val_length <- ifelse(is.vector(obj$standata[[param_mx]]), 1, ncol(obj$standata[[param_mx]]))
    } else {
        # If matrix not set, then val can be any length
        val_length <- length(val)
    }
    if (length(val)!=val_length){
        if (length(val)==1) {
            val <- rep(val, val_length)
        } else {
            if (!is.na(param_mx)){
                stop(paste0("Incorrect length for 'val'. Should have length 1 or length equal to the number of columns of ", param_mx))    
            } else {
                stop(paste0("Incorrect length for 'val'.  Should have length 1."))
            }
        }
    }
    obj$structure[[param]] <- val
    obj
}

##' @rdname outsim_update_parameter
##' @param df Degrees of freedom for exposure spline.
##' @param fn Name of spline-generating function. Defaults to \code{\link[splines]{ns}}
##' @param ... Additional arguments passed to \code{\link{create_spline_matrix}}.
##' @seealso \code{\link{create_spline_matrix}}
##' @export
outsim_update_exp_splines <- function(obj, df=1, fn="ns", ...){
    if (!inherits(obj, c("outsim"))) stop("'obj' must be of class 'outsim'.")
    if (all(is.na(obj$structure$x))) stop("obj$structure$x must exist.")
    Hx <- create_spline_matrix(obj$structure$x, df=df, fn=fn, ...)
    obj$standata <- add_Hx_standata(obj$standata, Hx)
    obj
}


expit <- function(x) {1/(1 + exp(-x))}
logit <- function(x) {log(x/(1-x))}