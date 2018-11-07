######################
# Contains the following functions
#
# expsim_create_study_skeleton_parallel()
# expsim_create_study_skeleton()
# expsim_update_times()
# expsim_update_time_splines()
# expsim_update_time_effect()
# expsim_update_parameter()
######################


##' @title Create List Structure for Simulations
##' @description These functions construct lists containing the structure for simulating exposure data and the objects needed for fitting the model to the data.
##' @param design String providing the study design. Currently only \code{"parallel"} is implemented.
##' @param ... Additional arguments passed to the design-specific functions.
##' @details These functions create the blank structure of groups, clusters, households, and subject.
##' Once parameters are set via \code{\link{expsim_update_parameter}} and related functions, then data can be sampled via
##' \code{\link{expsim_sample_observations}}
##' and posterior parameter estimates obtained from \code{\link{sample_exposure_model}}.
##' @return A list containing two sublists: \code{structure}, which contains settings and parameters for generating data, and \code{standata} which contains the study data as a \code{standata_exposure} obejct for sampling via \code{\link{sample_exposure_model}}.
#' @family exposure simulation functions
#' @seealso \code{\link{create_outcome_simulation_skeleton}}
#' @export
create_exposure_simulation_skeleton <- function(design="parallel",...){
    if (design=="parallel"){
        obj <- create_exposure_simulation_skeleton_parallel(...)
    } else if (design=="crossover"){
        stop("Crossover design not yet implemented.")
        # obj <- expsim_create_study_skeleton_crossover(...)
    } else {
        stop(paste0("Design = ", design, " is not supported. Please use 'parallel'."))
    }
    obj
}


##' @rdname create_exposure_simulation_skeleton
##' @param ngroups Number of groups.
##' @param nclusters Number of clusters. Either a single value or a vector of length \code{ngroups}.
##' @param nhouseholds Number of households in each cluster. If length 1, repeated for all clusters.
##' @param nobs Number of observations in each household. If length 1, repeated for all households and if length \code{sum(nclusters)}, repeated for all households within each cluster.
##' @param verbose Logical indicator for message printing.
##' @export
create_exposure_simulation_skeleton_parallel <- function(ngroups=1,
                                                  nclusters=1,
                                                  nhouseholds=1,
                                                  nobs=1,
                                                  verbose=TRUE){

    # Check input
    if(!is.null(nclusters) && !check_all_nonnegative_value(nclusters)) stop("'nclusters' should be a non-negative integer and must have at least one positive element.")
    if(!check_all_nonnegative_value(nhouseholds)) stop("'nhouseholds' should be a non-negative integer and must have at least one positive element.")
    if(!check_all_postive_value(ngroups)) stop("'ngroups' should be a positive integer.")

    # Process nclusters
    if(!length(nclusters) %in% c(1, ngroups)) stop("'nclusters' must have length 1 or 'ngroups'.")
    if (length(nclusters)==1) {
        nclusters <- rep(nclusters, times=ngroups)
        if(verbose) {
            message("Only one value of 'nclusters' provided. Applying value to all groups.")
        }
    }
    K <- sum(nclusters)

    # Process nhouseholds
    if(!length(nhouseholds) %in% c(1, K)) stop("'nhouseholds' must have length 1 or 'sum(nclusters)'.")
    if (length(nhouseholds)==1) {
        nhouseholds <- rep(nhouseholds, times=K)
        if(verbose) {
            message("Only one value of 'nhouseholds' provided. Applying value to all clusters.")
        }
    }
    H <- sum(nhouseholds)

    # Process nobs
    if(!length(nobs) %in% c(1, K, H)) stop("'nobs' must have length 1, 'sum(nclusters)', or 'sum(nhouseholds)'.")
    if (length(nobs)==1) {
        nobs <- rep(nobs, times=H)
        if(verbose) {
            message("Only one value of 'nobs' provided. Applying value to all households.")
        }
    } else if (length(nobs)==K) {
        nobs <- rep(nobs, times=nhouseholds)
        if(verbose) {
            message("One value of 'nobs' provided for each cluster. Applying each value to all households within corresponding cluster.")
        }
    }
    N <- sum(nobs)

    # Defaults
    group_of_cluster <- rep(1:ngroups, times=nclusters)
    cluster_of_hh <- rep(1:K, times=nhouseholds)
    hh_of_obs <- rep(1:H, times=nobs)
    cluster_of_obs <- cluster_of_hh[hh_of_obs]
    group_of_obs <- group_of_cluster[cluster_of_obs]

    study_structure <- list(ngroups=ngroups,
                            nclusters=nclusters,
                            nhouseholds=nhouseholds,
                            nobs=nobs,
                            group_of_cluster=group_of_cluster,
                            cluster_of_hh=cluster_of_hh,
                            hh_of_obs=hh_of_obs,
                            cluster_of_obs=cluster_of_obs,
                            group_of_obs=group_of_obs,
                            etaG=rep(NA, ngroups), # Group Means
                            sigK=rep(NA, ngroups), # Group-level standard deviations for drawing reK)
                            reK=rep(NA, K), # Cluster random effects
                            sigH=NA, # SD for household-level RE's. Not within cluster
                            reH=rep(NA, H), # Household-level random effects
                            sigW=NA, # SD for instrument/residual variation
                            times=NA,
                            timefn=function(t) {rep(0, length(t))},
                            meanW=rep(NA, N), # mean for generating data. will also have time component added.
                            design="parallel")
    study_standata <- list(G=ngroups,
                           K=K,
                           H=H,
                           N=N,
                           Ht=matrix(0, 0, 0),
                           times=NA,
                           timedf=0,
                           w=rep(NA, N),
                           group_of_cluster=group_of_cluster,
                           group_of_obs=group_of_obs,
                           cluster_of_obs=cluster_of_obs,
                           cluster_of_hh=cluster_of_hh,
                           hh_of_obs=hh_of_obs)

    class(study_standata) <- "standata_exposure"

    obj <- list(structure=study_structure,
         standata=study_standata)
    class(obj) <- "expsim"
    obj
}




#' @title Update Exposure Model Parameters
#' @description Sets or overwrites data-generating parameters for an exposure model object
#' @param obj 'expsim' object to update
#' @param param character string giving name of parameter to update. If not provided, parameter name is taken from \code{level} and \code{type}
#' @param level character string of model level being updated. See details.
#' @param type character string giving type of parameter being updated. See Details.
#' @param value value(s) to which parameter should be set
#' @param draw if \code{TRUE}, then the parameter value is sampled instead of being set to \code{val}
#' @details Only specific combinations of \code{level} and \code{type} are allowed. They are:
#' \itemize{
#' \item \code{"group"}: \code{"mean"} sets 'etaG'
#' \item \code{"cluster"}: \code{"sd"} and \code{"re"} set 'sigK' and 'reK', respectively
#' \item \code{"household"}: \code{"sd"} and \code{"re"} set 'sigH' and 'reH', respectively
#' \item \code{"instrument"}: \code{"sd"} sets 'sigW'.
#' \item \code{"time"}: \code{"mean"} sets 'timefn'.
#' }
#' @family exposure simulation functions
#' @export
#' @importFrom stats rnorm
expsim_update_parameter <- function(obj, param, level=c("group", "cluster","household", "instrument", "correlation","time"), type=c("mean", "sd", "re"), value=NULL, draw=is.null(value)){

    if (!inherits(obj, "expsim")) stop("'obj' must be of class 'expsim'.")
    if (!missing(param)){
        param_list <- c("etaG", "sigK","reK", "sigH", "reH",  "sigW", "timefn")
        param_check <- param %in% param_list
        if (!param_check) stop('"param" must be one of c("etaG","sigK", "reK", "sigH", "reH",  "sigW", "timefn")')
    } else {
        level <- match.arg(level)
        type <- match.arg(type)
        param <- switch(paste0(level, "_", type),
               group_mean="etaG",
               group_sd=NA,
               group_re=NA,
               cluster_mean=NA,
               cluster_sd="sigK",
               cluster_re="reK",
               household_mean=NA,
               household_sd="sigH",
               household_re="reH",
               instrument_mean=NA,
               instrument_sd="sigW",
               instrument_re=NA,
               time_mean="timefn",
               time_sd=NA,
               time_re=NA,
               NA)
        if (is.na(param)) stop(paste0(level, "_", type, " is not a supported parameter."))
    }

    # Need to do some error/length checking.
    current_length <- length(obj$structure[[param]])
    if (draw){
        if(!is.null(value)) warning("'value' provided but 'draw=TRUE'. Ignoring provided value of 'value'.")

        if (!param %in% c("reK", "reH")) {
            stop(paste0(param, " should be set, not drawn from a distribution."))
        }

        # Need to check that mean/sd has been set, first.
        cur_mean <- switch(param,
                           reK=0,
                           reH=0)
        cur_sd <- switch(param,
                         reK=obj$structure$sigK,
                         reH=obj$structure$sigH)
        cur_n <- switch(param,
                        reK=obj$standata$K,
                        reH=obj$standata$H)
        value <- rnorm(cur_n, cur_mean, cur_sd)
    }
    if (length(value)!=current_length) {
        if (length(value)==1){
            value <- rep(value, current_length)
        } else {
            stop(paste0("'value' should have length ", current_length, "."))
        }
    }
    obj$structure[[param]] <- value

    obj
}

##' @title Update time values and spline object
##' @description Creates times and splines of time for simulations
#' @param obj Simualtion object, created by \code{\link{create_exposure_simulation_skeleton}} or \code{\link{create_outcome_simulation_skeleton}}.
##' @param df degrees of freedom for the time spline used in model estimation.
##' @param fn function for generating splines. Defaults to \code{\link{ns}}
##' @param ... Additional arguments passed to \code{}
#' @family exposure simulation functions
#' @family outcome simulation functions
##' @seealso \code{\link{create_spline_range}} \code{\link{add_Ht_standata}}
##' @export
##' @importFrom splines ns
sim_update_time_splines <- function(obj, df=1, fn="ns", ...){
    if (!inherits(obj, c("expsim", "outsim"))) stop("'obj' must be of class 'expsim' or 'outsim'.")
    if (all(is.na(obj$structure$times))) stop("obj$structure$times must exist.")
    Ht <- create_spline_matrix(obj$structure$times, df=df, fn=fn, ...)
    obj$standata <- add_Ht_standata(obj$standata, Ht)
    obj
}


#' @title Sample Simulated Exposure Observations
#' @description Draws a sample of exposure observations using the specified model parameters
#' @param obj exposure simualtion object, created by \code{\link{create_exposure_simulation_skeleton}}.
#' @details This function assumes that \code{etaK} (or at least \code{etaG}), \code{sigW}, \code{reH}, and \code{timefn} have been set. The latter defaults to zero, but the others must be set with \code{\link{expsim_update_parameter}}. If these have not been set, then this will set \code{w} to a vector of \code{NA}.
#' @family exposure simulation functions
#' @export
#' @importFrom stats rnorm
expsim_sample_observations <- function(obj){
        if (!inherits(obj, "expsim")) stop("'obj' must be of class 'expsim'.")

    obj$structure$meanW <- with(obj$structure, etaG[group_of_obs] + reH[hh_of_obs] + timefn(times))

    if (!any(is.na(obj$structure$reK))){
        obj$structure$meanW <- obj$structure$meanW + obj$structure$reK[obj$structure$cluster_of_obs]
    }
    obj$standata$w=rnorm(n=obj$standata$N, mean=obj$structure$meanW, sd=obj$structure$sigW)
    obj
}




## Essentially same as other function, but kept separate since runif()
## Perhaps put into an 'update_data' function in future.
##' @rdname sim_update_time_splines
##' @param times Values to set as the times.
##' @param draw Logical indicating times should be sampled from a Unif(0,1) distribution.
##' @export
##' @importFrom stats runif
sim_update_times <- function(obj, times=NULL, draw=is.null(times)){
        if (!inherits(obj, c("expsim", "outsim"))) stop("'obj' must be of class 'expsim' or 'outsim'.")

    if (draw){
        if(!is.null(times)) warning("'times' is provided, but 'draw=TRUE'. Ignoring provided value of 'times'.")
        # Draw from (0, 1)
        times <- stats::runif(obj$standata$N)
    }
    if (is.null(times)) {
        stop("'times' is NULL, but 'draw=FALSE'.")
    }
    obj$structure$times <- times
    obj$standata$times <- times

    obj
}




