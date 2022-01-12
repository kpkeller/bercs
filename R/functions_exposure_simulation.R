######################
# Contains the following functions
#
# create_exposure_simulation_skeleton()
# create_exposure_simulation_skeleton_parallel()
# expsim_update_times()
# sim_update_time_splines()
# expsim_update_time_effect()
# expsim_update_parameter()
# expsim_sample_observations()
# sim_update_times()
######################


##' @title Create List Structure for Simulations
##' @description These functions construct lists containing the structure for simulating exposure data and the objects needed for fitting the model to the data.
##' @param design String providing the study design. Currently only \code{"parallel"} and \code{"crossover"} is implemented.
##' @param ... Additional arguments passed to the design-specific functions.
##' @details These functions create the blank structure of groups, clusters, and units (e.g. households).
##' Once parameters are set via \code{\link{expsim_update_parameter}} and related functions, then data can be sampled via
##' \code{\link{expsim_sample_observations}}
##' and posterior parameter estimates obtained from \code{\link{sample_exposure_model}}.
##' See \code{\link{sample_exposure_model}} for a mathematical description of the model.
##' @return A list containing two sublists: \code{structure}, which contains settings and parameters for generating data, and \code{standata} which contains the study data as a \code{standata_exposure} obejct for sampling via \code{\link{sample_exposure_model}}.
#' @family exposure simulation functions
#' @seealso \code{\link{create_outcome_simulation_skeleton}}, \code{\link{sample_exposure_model}}
#' @export
create_exposure_simulation_skeleton <- function(design="parallel",...){
    if (design=="parallel"){
        obj <- create_exposure_simulation_skeleton_parallel(...)
    } else if (design=="crossover"){
        obj <- create_exposure_simulation_skeleton_crossover(...)
    } else {
        stop(paste0("Design = ", design, " is not supported. Please use 'parallel' or 'crossover'."))
    }
    obj
}


##' @rdname create_exposure_simulation_skeleton
##' @param ngroups Number of groups.
##' @param nclusters Number of clusters. Either a single value or a vector of length \code{ngroups}.
##' @param nunits Number of units (e.g. households) in each cluster. If length 1, repeated for all clusters.
##' @param nobs Number of observations in each unit. If length 1, repeated for all units and if length \code{sum(nclusters)}, repeated for all units within each cluster.
##' @param etaG Group means. Optional; can be set later with \code{\link{expsim_update_parameter}}.
##' @param sigI Standard deviation for household random effect. Optional; can be set later with \code{\link{expsim_update_parameter}}.
##' @param sigW Standard deviation of residual variation for each observation. Optional; can be set later with \code{\link{expsim_update_parameter}}.
##' @param verbose Logical indicator for message printing.
##' @param time optional times for observations. Can be set later via \code{\link{sim_update_times}}.
##' @param timefn Function relating times to mean exposure. Can be set later via \code{\link{expsim_update_parameter}}.
##' @examples
##' # Create simulation structure for parallel design
##' es <- create_exposure_simulation_skeleton_parallel(ngroups=2, nunits=3)
##' str(es)
##' @export
create_exposure_simulation_skeleton_parallel <- function(ngroups=1,
                                                  nclusters=1,
                                                  nunits=1,
                                                  nobs=1,
                                                  etaG=rep(NA, ngroups),
                                                  sigI=NA,
                                                  sigW=NA,
                                                  time=NA,
                                                  timefn=function(t) {rep(0, length(t))},
                                                  verbose=TRUE){

    # Check input
    if(!is.null(nclusters) && !check_all_nonnegative_value(nclusters)) stop("'nclusters' should be a non-negative integer and must have at least one positive element.")
    if(!check_all_nonnegative_value(nunits)) stop("'nunits' should be a non-negative integer and must have at least one positive element.")
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

    # Process nunits
    if(!length(nunits) %in% c(1, K)) stop("'nunits' must have length 1 or 'sum(nclusters)'.")
    if (length(nunits)==1) {
        nunits <- rep(nunits, times=K)
        if(verbose) {
            message("Only one value of 'nunits' provided. Applying value to all clusters.")
        }
    }
    n <- sum(nunits)

    # Process nobs
    if(!length(nobs) %in% c(1, K, n)) stop("'nobs' must have length 1, 'sum(nclusters)', or 'sum(nunits)'.")
    if (length(nobs)==1) {
        nobs <- rep(nobs, times=n)
        if(verbose) {
            message("Only one value of 'nobs' provided. Applying value to all units")
        }
    } else if (length(nobs)==K) {
        nobs <- rep(nobs, times=nunits)
        if(verbose) {
            message("One value of 'nobs' provided for each cluster. Applying each value to all units within corresponding cluster.")
        }
    }
    N <- sum(nobs)

    # Defaults
    group_of_cluster <- rep(1:ngroups, times=nclusters)
    cluster_of_unit <- rep(1:K, times=nunits)
    unit_of_obs <- rep(1:n, times=nobs)
    cluster_of_obs <- cluster_of_unit[unit_of_obs]
    group_of_obs <- group_of_cluster[cluster_of_obs]

    study_structure <- list(ngroups=ngroups,
                            nclusters=nclusters,
                            nunits=nunits,
                            nobs=nobs,
                            group_of_cluster=group_of_cluster,
                            cluster_of_unit=cluster_of_unit,
                            unit_of_obs=unit_of_obs,
                            cluster_of_obs=cluster_of_obs,
                            group_of_obs=group_of_obs,
                            etaG=etaG, # Group Means
                            sigK=rep(NA, ngroups), # Group-level standard deviations for drawing reK)
                            reK=rep(NA, K), # Cluster random effects
                            sigI=sigI, # SD for unit-level RE's. Not within cluster
                            reI=rep(NA, n), # unit-level random effects
                            sigW=sigW, # SD for observation residual variation
                            time=time,
                            timefn=timefn,
                            meanW=rep(NA, N), # mean for generating data. will also have time component added.
                            design="parallel")
    study_standata <- list(G=ngroups,
                           K=K,
                           n=n,
                           N=N,
                           Mt=matrix(0, 0, 0),
                           time=NA,
                           timedf=0,
                           w=rep(NA, N),
                           group_of_cluster=group_of_cluster,
                           group_of_obs=group_of_obs,
                           cluster_of_obs=cluster_of_obs,
                           cluster_of_unit=cluster_of_unit,
                           unit_of_obs=unit_of_obs)

    class(study_standata) <- "standata_exposure"

    obj <- list(structure=study_structure,
         standata=study_standata)
    class(obj) <- "expsim"
    obj
}


##' @rdname create_exposure_simulation_skeleton
##' @param nobs Number of observations in each unit. If length 1, repeated for all units and if length \code{sum(nclusters)}, repeated for all units within each cluster.
##' @param nobs1 Number of observations in the first group for each unit (should have length equal to the number of units). Note that the "first" group might not be group 1.
##' @param nobs2 Number of observations in the second group for each unit (should have length equal to the number of units). Note that the "second" group might not be group 2.
##' @param firstgroup Integer indicating which group a unit is in first (should have length equal to the number of units).
##' @details The crossover design is currently implemented with only 2 groups.
##' @examples
##' # Create simulation structure for crossover design
##' # with 10 units and 3 measurements per unit under
##' # each treatment group (60 measures total)
##' es <- create_exposure_simulation_skeleton_crossover(ngroups=2, nunits=10, nobs1=3, nobs2=3, firstgroup=rep(1:2, each=5))
##' str(es)
##' @export
create_exposure_simulation_skeleton_crossover <- function(ngroups=2,
                                                         nclusters=1,
                                                         nunits=1,
                                                         nobs1=1,
                                                         nobs2=1,
                                                         firstgroup=1,
                                                         etaG=rep(NA, ngroups),
                                                         sigI=NA,
                                                         sigW=NA,
                                                         time=NA,
                                                         timefn=function(t) {rep(0, length(t))},
                                                         verbose=TRUE){

    # Check input
    if(!is.null(nclusters) && !check_all_nonnegative_value(nclusters)) stop("'nclusters' should be a non-negative integer and must have at least one positive element.")
    if(!check_all_nonnegative_value(nunits)) stop("'nunits' should be a non-negative integer and must have at least one positive element.")
    if(!check_all_postive_value(ngroups)) stop("'ngroups' should be a positive integer.")

    if (ngroups!=2) stop("Currently only 'ngroups=2' is supported for crossover designs.")

    # Process nclusters
    if (nclusters!=1) stop("Multiple clusters not currently supported for crossover design.")
    # if(!length(nclusters) %in% c(1, ngroups)) stop("'nclusters' must have length 1 or 'ngroups'.")
    # if (length(nclusters)==1) {
    #     nclusters <- rep(nclusters, times=ngroups)
    #     if(verbose) {
    #         message("Only one value of 'nclusters' provided. Applying value to all groups.")
    #     }
    # }
    K <- sum(nclusters)

    # Process nunits
    # if(!length(nunits) %in% c(1, K)) stop("'nunits' must have length 1 or 'sum(nclusters)'.")
    # if (length(nunits)==1) {
    #     nunits <- rep(nunits, times=K)
    #     if(verbose) {
    #         message("Only one value of 'nunits' provided. Applying value to all clusters.")
    #     }
    # }
    n <- sum(nunits)

    # Process nobs
    if(!length(nobs1) %in% c(1, n)) stop("'nobs1' must have length 1 or 'sum(nunits)'.")
    if(!length(nobs2) %in% c(1, n)) stop("'nobs2' must have length 1 or 'sum(nunits)'.")
    if (length(nobs1)==1) {
        nobs1 <- rep(nobs1, times=n)
        if(verbose) {
            message("Only one value of 'nobs1' provided. Applying value to all units")
        }
    }
    if (length(nobs2)==1) {
        nobs2 <- rep(nobs2, times=n)
        if(verbose) {
            message("Only one value of 'nobs2' provided. Applying value to all units")
        }
    }
    N <- sum(nobs1) + sum(nobs2)


    # Process firstgroup
    if(!length(firstgroup) %in% c(1, n)) stop("'firstgroup' must have length 1 or 'sum(nunits)'.")
    if (length(firstgroup)==1) {
        firstgroup <- rep(firstgroup, times=n)
        if(verbose) {
            message("Only one value of 'firstgroup' provided. Applying value to all units")
        }
    }



    # Defaults
    # group_of_cluster <- rep(1:ngroups, times=nclusters)
    group_of_cluster <- NA
    cluster_of_unit <- rep(1, times=n) #rep(1:K, times=nunits)
    unit_of_obs <- rep(1:n, times=nobs1 + nobs2)
    cluster_of_obs <- cluster_of_unit[unit_of_obs]
    group_of_obs <- rep(c(rbind(firstgroup,
                                3- firstgroup)),
                        times=c(rbind(nobs1,
                                      nobs2))) #group_of_cluster[cluster_of_obs]

    study_structure <- list(ngroups=ngroups,
                            nclusters=nclusters,
                            nunits=nunits,
                            nobs1=nobs1,
                            nobs2=nobs2,
                            nobs=nobs1 + nobs2,
                            group_of_cluster=group_of_cluster,
                            cluster_of_unit=cluster_of_unit,
                            unit_of_obs=unit_of_obs,
                            cluster_of_obs=cluster_of_obs,
                            group_of_obs=group_of_obs,
                            etaG=etaG, # Group Means
                            sigK=rep(NA, ngroups), # Group-level standard deviations for drawing reK)
                            reK=rep(NA, K), # Cluster random effects
                            sigI=sigI, # SD for unit-level RE's. Not within cluster
                            reI=rep(NA, n), # unit-level random effects
                            sigW=sigW, # SD for observation residual variation
                            time=time,
                            timefn=timefn,
                            meanW=rep(NA, N), # mean for generating data. will also have time component added.
                            design="crossover")
    study_standata <- list(G=ngroups,
                           K=K,
                           n=n,
                           N=N,
                           Mt=matrix(0, 0, 0),
                           time=NA,
                           timedf=0,
                           w=rep(NA, N),
                           group_of_cluster=group_of_cluster,
                           group_of_obs=group_of_obs,
                           cluster_of_obs=cluster_of_obs,
                           cluster_of_unit=cluster_of_unit,
                           unit_of_obs=unit_of_obs)

    class(study_standata) <- "standata_exposure"

    obj <- list(structure=study_structure,
                standata=study_standata)
    class(obj) <- "expsim"
    obj
}







#' @title Update Exposure Model Parameters
#' @description Sets or overwrites data-generating parameters for an exposure model object
#' @param obj 'expsim' object to update (see \code{\link{create_exposure_simulation_skeleton}}).
#' @param param character string giving name of parameter to update. If not provided, parameter name is taken from \code{level} and \code{type}.
#' @param level character string of model level being updated. Must be one of "group", "cluster", "unit", "observation", or "time". See Details.
#' @param type character string giving type of parameter being updated. See Details.
#' @param value value(s) to which parameter should be set
#' @param draw if \code{TRUE}, then the parameter value is sampled instead of being set to \code{value}
#' @details Only specific combinations of \code{level} and \code{type} are allowed. They are:
#' \itemize{
#' \item \code{"group"}: \code{"mean"} sets 'etaG'
#' \item \code{"cluster"}: \code{"sd"} and \code{"re"} set 'sigK' and 'reK', respectively
#' \item \code{"unit"}: \code{"sd"} and \code{"re"} set 'sigI' and 'reI', respectively
#' \item \code{"observation"}: \code{"sd"} sets 'sigW'.
#' \item \code{"time"}: \code{"fn"} sets 'timefn'. To update the times themselves, use \code{\link{sim_update_times}}.
#' }
#' After setting and/or sampling values via this funciton, use \code{\link{expsim_sample_observations}} to sample the exposure values.
#' @family exposure simulation functions
#' @export
#' @importFrom stats rnorm
#' @examples
#' # Create simulation structure
#' es <- create_exposure_simulation_skeleton_parallel(ngroups=2, nunits=10)
#'
#' # No standard deviation parameter set by default
#' es$structure$sigK
#' # Add standard deviation parameter value
#' es <- expsim_update_parameter(es, level="cluster", type="sd", value=1)
#' es$structure$sigK
#'
#' # Sample random effect values
#' es <- expsim_update_parameter(es, level="cluster", type="re", draw=TRUE)
#'
#' # Add time function
#' es <- expsim_update_parameter(es, level="time", type="mean", value=function(w) 2*cos(w))
expsim_update_parameter <- function(obj, param, level=c("group", "cluster","unit", "observation", "correlation","time"), type=c("mean", "sd", "re", "fn"), value=NULL, draw=is.null(value)){

    if (!inherits(obj, "expsim")) stop("'obj' must be of class 'expsim'.")
    if (!missing(param)){
        param_list <- c("etaG", "sigK","reK", "sigI", "reI",  "sigW", "timefn")
        param_check <- param %in% param_list
        if (!param_check) stop('"param" must be one of c("etaG","sigK", "reK", "sigI", "reI",  "sigW", "timefn")')
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
               unit_mean=NA,
               unit_sd="sigI",
               unit_re="reI",
               observation_mean=NA,
               observation_sd="sigW",
               observation_re=NA,
               time_fn="timefn",
               time_sd=NA,
               time_re=NA,
               NA)
        if (is.na(param)) stop(paste0(level, "_", type, " is not a supported parameter."))
    }

    # Need to do some error/length checking.
    current_length <- length(obj$structure[[param]])
    if (draw){
        if(!is.null(value)) warning("'value' provided but 'draw=TRUE'. Ignoring provided value of 'value'.")

        if (!param %in% c("reK", "reI")) {
            stop(paste0(param, " should be set, not drawn from a distribution."))
        }

        # Need to check that mean/sd has been set, first.
        cur_mean <- switch(param,
                           reK=0,
                           reI=0)
        cur_sd <- switch(param,
                         reK=obj$structure$sigK,
                         reI=obj$structure$sigI)
        cur_n <- switch(param,
                        reK=obj$standata$K,
                        reI=obj$standata$n)
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
#' @param obj Simulation object, created by \code{\link{create_exposure_simulation_skeleton}} or \code{\link{create_outcome_simulation_skeleton}}.
##' @param time Values to set as the times. Must have length equal to the number of observations in \code{obj}.
##' @param draw Logical indicating times should be sampled from a Unif(0,1) distribution.
##' @param ... Passed to \code{runif}
##' @export
##' @importFrom stats runif
##' @family exposure simulation functions
##' @family outcome simulation functions
##' @seealso \code{\link{create_spline}} \code{\link{add_spline_time}}
##' @examples
##' # Create structure
##' es <- create_exposure_simulation_skeleton_parallel(ngroups=2, nunits=10)
##'
##' es <- sim_update_times(es, draw=TRUE)
##' es <- sim_update_time_splines(es, df=3)
sim_update_times <- function(obj, time=NULL, draw=is.null(time), ...){
    if (!inherits(obj, c("expsim", "outsim"))) stop("'obj' must be of class 'expsim' or 'outsim'.")

    if (draw){
        if(!is.null(time)) warning("'time' is provided, but 'draw=TRUE'. Ignoring provided value of 'time'.")
        # Draw from (0, 1)
        time <- stats::runif(obj$standata$N,...)
    }
    if (is.null(time)) {
        stop("'time' is NULL, but 'draw=FALSE'.")
    }
    if (length(time) !=obj$standata$N){
        stop("Incorrect length for times.")
    }
    obj$structure$time <- time
    obj$standata$time <- time

    obj
}

##' @rdname sim_update_times
##' @param df degrees of freedom for the time spline used in model estimation.
##' @param fn function for generating splines. Defaults to \code{\link{ns}}.
##' @param ... Additional arguments passed to \code{fn()} via \code{\link{create_spline_matrix}}.
##' @export
##' @importFrom splines ns
sim_update_time_splines <- function(obj, df=1, fn="ns", ...){
    if (!inherits(obj, c("expsim", "outsim"))) stop("'obj' must be of class 'expsim' or 'outsim'.")
    if (all(is.na(obj$structure$time))) stop("obj$structure$time must exist.")
    Mt <- create_spline_matrix(obj$structure$time, df=df, fn=fn, ...)
    obj$standata <- add_spline_time(obj$standata, Mt)
    obj
}


#' @title Sample Simulated Exposure Observations
#' @description Draws a sample of exposure observations using the specified model parameters
#' @param obj exposure simulation object, created by \code{\link{create_exposure_simulation_skeleton}}.
##' @details This function assumes that \code{etaG}, \code{sigW}, \code{reI}, and \code{timefn} have been set. Optionally, \code{sigK} should be set if a cluster random effect is to be included. These can be set with \code{\link{expsim_update_parameter}}. Using these values, it first computes the (conditional) observation mean, with is
##' @return An object of class \code{expsim}, with updated values of \code{obj$standata$w} and \code{obj$structure$meanW}.
##' @family exposure simulation functions
##' @export
##' @importFrom stats rnorm
##' @examples
##' es <- create_exposure_simulation_skeleton_parallel(ngroups=2,
##' nclusters=15,
##' nunits=10,
##' nobs=2)
##' # Set group means
##' es <- expsim_update_parameter(es,
##'                               level="group",
##'                               type="mean",
##'                               value=c(3, 4))
##' # Set unit random effect standard deviation
##' es <- expsim_update_parameter(es,
##'                               level="unit",
##'                               type="sd",
##'                               value=c(1))
##' # Set observation standard deviation
##' es <- expsim_update_parameter(es,
##'                               level="observation",
##'                               type="sd",
##'                               value=c(1))
##' # Sample unit random effects
##' es <- expsim_update_parameter(es,
##'                               level="unit",
##'                               type="re")
##' # Compute (conditional) observation means and sample observations
##' expsim_demo <- expsim_sample_observations(es)
##' str(expsim_demo$structure$meanW)
##' str(expsim_demo$standata$w)
##' # Add priors and estimate parameters
##' expsim_demo$standata <- add_priors(expsim_demo$standata)
##' fit <- sample_exposure_model(expsim_demo$standata)
##' print(fit, pars=c("reI_raw", "muW", "reI", "reK", "reK_raw"), include=FALSE)
expsim_sample_observations <- function(obj){
    if (!inherits(obj, "expsim")) stop("'obj' must be of class 'expsim'.")
    if (any(is.na(obj$structure$etaG))) stop("etaG must be set. See expsim_update_parameter().")
    if (any(is.na(obj$structure$reI))) stop("reI must be set. See expsim_update_parameter().")
    if (any(is.null(obj$structure$timefn))) stop("timefn must be set. See expsim_update_parameter().")
    if (any(is.na(obj$structure$sigW))) stop("sigW must be set. See expsim_update_parameter().")
    obj$structure$meanW <- with(obj$structure, etaG[group_of_obs] + reI[unit_of_obs] + timefn(time))

    if (!any(is.na(obj$structure$reK))){
        obj$structure$meanW <- obj$structure$meanW + obj$structure$reK[obj$structure$cluster_of_obs]
    }
    obj$standata$w=rnorm(n=obj$standata$N, mean=obj$structure$meanW, sd=obj$structure$sigW)
    obj
}





