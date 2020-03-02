# Functions to support fitting exposure models
#
# create_standata_exposure()
# sample_exposure_model()
# compute_fitted_mean()
# plot_exposure_means_bytime()
# plot_exposure_means_boxplot()

##' @title Create List for Fitting Exposure Model via STAN
##' @description Creates a list object that follows the structure required by
##'  the exposure models that are implemented in STAN.
##' @param data Optional data frame containing data in a long format.
##' @param group Vector of group assignments. See Details.
##' @param conc Vector of exposure concentrations.
##' @param unit_id Vector identifying the distinct unit (e.g. person, household) for each observation.
##' @param clust_id Vector identifying cluster membership for each observation.
##' @param time Vector of times corresponding to \code{conc} values.
##' @param Mt Matrix of splines values for time to include in the model. Defaults to NULL, and can be added later via \code{\link{add_spline_time}}.
##' @param log_transform Should the concentration value be log-transformed?
##' @param return_addition See 'Value'.
##' @details This function takes as input a `long` data frame and extracts from it the information needed for fitting the exposure model in STAN.
##'
##' The values of \code{group}, \code{unit_id}, and \code{clust_id} are converted to factor and then coerced to integer. The ordering of groups (units, clusters) in the output object will depend on the default ordering introduced by \code{\link{factor}}.
##' @return A list of class \code{standata_exposure} that contains the following:
##' \itemize{
##' \item {\code{G} -- Number of groups}
##' \item {\code{K} -- Number of clusters.}
##' \item{\code{H} Number of households}
##' \item{\code{N} -- Number of observations}
##' \item{\code{cluster_of_obs} -- Integer providing cluster number of each observation}
##' \item{\code{group_of_obs} -- Integer providing group number of each observation}
##' \item{\code{hh_of_obs} -- Integer providing household number of each observation}
##' \item{\code{w} -- Concentration observations}
##' \item{\code{times} -- The date/time of each observation. This is not used directly in model fitting, but contained in the object to facilitate plotting and other summaries.}
##' }
##' If \code{return_addition=TRUE}, then a two-element list is returned. The first element
##' is the \code{standata_exposure} object described above. The second element is a modified version of \code{data}, with the variables \code{group_of_obs}, \code{cluster_of_obs}, \code{hh_of_obs}, and \code{times} added (or overwritten).
##' @seealso \code{\link{sample_exposure_model}}, \code{\link{create_standata_outcome}}
##' @export
create_standata_exposure <- function(data=NULL,
                                     group=data$group,
                                     conc=data$conc,
                                     unit_id=data$unit_id,
                                     clust_id=data$clust_id,
                                     time=data$time,
                                     Mt=NULL,
                                     log_transform=FALSE,
                                     return_addition=FALSE){
    N <- length(conc)
    if(!all(length(group)==N,
           length(unit_id)==N)) stop("The lengths of 'group', 'unit_id', and 'conc' must be equal")

    # Create data frame
    newdata <- data.frame(group_of_obs=as.numeric(factor(group)),
                          unit_of_obs=as.numeric(factor(unit_id)),
                          conc=conc)
    if (!is.null(clust_id) && any(clust_id!=0)){
        if (length(clust_id)!=N) stop("The lengths of 'clust_id' and 'conc' must be equal")
        newdata$cluster_of_obs <- as.numeric(factor(clust_id))
    } else {
        newdata$cluster_of_obs <- rep(0, N)
    }
    if (!is.null(time)){
        if (length(time)!=N) stop("The lengths of 'clust_id' and 'conc' must be equal")
        newdata$time <- time
    } else {
        newdata$time <- rep(0, N)
    }

    out <- list()

    out$G <- max(newdata$group_of_obs)
    out$K <- max(newdata$cluster_of_obs)
    out$H <- max(newdata$unit_of_obs)
    out$N <- N
    out$cluster_of_obs <- newdata$cluster_of_obs
    out$group_of_obs <-  newdata$group_of_obs
    out$unit_of_obs <- newdata$hh_of_obs

    if (log_transform){
        out$w <- log(newdata$conc)
    } else {
        out$w <- newdata$conc
    }
    out$time <- newdata$time
    out$timedf <- 0
    if (!is.null(Mt)){
        out <- add_spline_time(out, Ht)
    } else {
        out <- add_spline_time(out, array(0, dim=c(0, 0)))
    }
    class(out) <- "standata_exposure"
    if (return_addition){
        return(list(standata=out,
                    df=newdata))
    } else {
        return(out)
    }
}

##' @title Sample Exposure Model
##' @description Samples from the posterior distribution of the exposure model, using STAN
##' @param standata An object of class `standata_exposure`, typically created from \code{\link{create_standata_exposure}}.
##' @param B Number of post-warmup iterations.
##' @param warmup Number of warmup iterations.
##' @param chains Number of chains to sample.
##' @param control List provided as the \code{control} argument of \code{\link[rstan]{sampling}}
##' @param ... Additional arguments passed to \code{\link[rstan]{sampling}}.
##' @details The model is:
##' \deqn{w = \eta_G...}
##' @seealso \code{\link{create_standata_exposure}}, \code{\link{sample_outcome_model}}
##' @importFrom rstan sampling
##' @export
sample_exposure_model <- function(standata,
                                  B=1000,
                                  warmup=B,
                                  chains=4,
                                  control=list(adapt_delta=0.9,
                                               max_treedepth=12),
                                  ...){

    if (!inherits(standata, "standata_exposure")) stop("`standata` must be of class 'standata_exposure'.")

    model_name <- "exposure_model"
    exp_stanfit <- rstan::sampling(stanmodels[[model_name]],
                          data = standata,
                          iter = B + warmup,
                          warmup = warmup,
                          chains = chains,
                          control = control,
                          ...)
    exp_stanfit
}




#' @title Compute fitted concentration
#' @description Calculates the posterior mean exposure concentrations.
#' @param stanfit Fitted STAN model object containing posterior samples of parameters
#' @param standata List containing model structure information corresponding to \code{stanfit}.
#'  Typically a list of class \code{standata_exposure}. However, it may be any list (or data frame) containing
#'  \code{cluster_of_obs} (if \code{etaK} was sampled in \code{stanfit}), \code{group_of_obs}
#'  (if \code{etaK} was not sampled in \code{stanfit}),  \code{hh_of_obs} (if \code{include_reH=TRUE}),
#'   and \code{Ht} (if \code{include_time=TRUE}). This may be useful when calculating fitted means for records
#'    at different times from the observations used to fit the model.
#' @param parvalues List of parameter values to use for computing means. Not needed if \code{stanfit} is provided.
#' @param include_time Logical indicator of whether time should be included.
#' @param include_reH Logical indicator of whether the household-level random effect should be included.
#' @param exp_transform Logical indicator of whether concentrations should be exponentiated.
#' @param add Logical indicator of whether the modeled means should be added to \code{standata} or just returned directly (the default).
#' @param ... Additional arugments passed to \code{\link[rstan]{extract}}.
#' @seealso \code{\link{plot_exposure_means_bytime}}
#' @export
#' @importFrom rstan extract
compute_fitted_mean <- function(stanfit,
                                standata,
                                parvalues=NULL,
                                include_time=TRUE,
                                include_reH=TRUE,
                                exp_transform=FALSE,
                                add=FALSE,
                                ...){


    if (is.null(parvalues)){
        nocluster <- ifelse("reK[1]" %in% names(stanfit), FALSE, TRUE)

        pars <- "etaG"
        if (!nocluster) pars <- c(pars, "reK")
        if(include_reH) pars <- c(pars, "reH")
        if (include_time) pars <- c(pars, "theta")

        postmean <- rstan::extract(stanfit,
                             pars=pars, ...)
        postmean <- lapply(postmean, colMeans)
    } else {
        nocluster <- ifelse("reK" %in% names(parvalues), FALSE, TRUE)

        pars <- "etaG"
        if (!nocluster) pars <- c(pars, "reK")
        if(include_reH) pars <- c(pars, "reH")
        if (include_time) pars <- c(pars, "theta")
        postmean <- parvalues
        if (!all(pars %in% names(postmean))) stop("'parvalues' provided, but not all names present.")
    }

    ltmean <- postmean$etaG[standata$group_of_obs]
    if (!nocluster){
        ltmean <- ltmean + postmean$reK[standata$cluster_of_obs]
    }
    if (include_reH){
        ltmean <- ltmean + postmean$reH[standata$hh_of_obs]
    }
    if (include_time){
        ltmean <- ltmean + as.vector(standata$Ht %*% postmean$theta)
    }

    if (exp_transform) ltmean <- exp(ltmean)
    if (add) {
        standata$ltmean <- ltmean
        return(standata)
    } else {
        return(ltmean)
    }
}


##' @title Plot Exposure Model Summaries
##' @description Different plots of the concentrations modelled in the exposure model
#' @inheritParams compute_fitted_mean
#' @param group_names Names of groups to use in plot labels
#' @details Uses \code{\link[ggplot2]{ggplot}} to create plots of concentrations using posterior means of model parameters. The graphical object is returned and can be customized if needed.
##' @export
##' @import ggplot2
plot_exposure_means_bytime <- function(stanfit,
                                  standata,
                                  include_time=TRUE,
                                include_reH=TRUE,
                                  exp_transform=FALSE,
                                  group_names=paste0("Group ", 1:standata$G)){

    ltmean <-  compute_fitted_mean(stanfit = stanfit,
                                     standata = standata,
                                     include_time =include_time,
                                     include_reH=include_reH,
                                     add = FALSE,
                                     exp_transform = exp_transform)

    group <- as.factor(standata$group_of_obs)
    if(exp_transform) {
        wplot <-   exp(standata$w)
    } else {
        wplot <- standata$w
    }
    g <- ggplot()  + geom_point(aes(x=standata$times, y=wplot), col="black", shape=2)
    g <- g + geom_point(aes(x=standata$times, y=ltmean, col=group))
    g <- g + scale_color_discrete(name="Group", labels=group_names) + theme_bw() + xlab("Time") + ylab("Concentration")
    g
}

## NOTE-- this includes repeated measures, since it has one value per observation
##' @rdname plot_exposure_means_bytime
##' @param one_per_person Logical indicator of whether a single observation per person should be included. If \code{FALSE} (the default), then
##' all observations are included in the boxplot.
##' @export
##' @import ggplot2
plot_exposure_means_boxplot <- function(stanfit,
                                 standata,
                                 include_time=TRUE,
                                 include_reH=TRUE,
                                 exp_transform=FALSE,
                                 one_per_person=!include_time,
                                 group_names=paste0("Group ", 1:standata$G)){

   ltmean <-  compute_fitted_mean(stanfit = stanfit,standata = standata,
                                    include_time =include_time,
                                    include_reH=include_reH,
                                    add = FALSE,
                                    exp_transform = exp_transform)

    group <- as.factor(standata$group_of_obs)
   if (one_per_person){
       inds <- !duplicated(standata$hh_of_obs)
       ltmean <- ltmean[inds]
       group <- group[inds]
   }
    g <- ggplot() + geom_boxplot(aes(x=group, y=ltmean, fill=group))  + scale_fill_discrete(name="Group", labels=group_names) + theme_bw() + xlab("Time") + ylab("Concentration")
    g
}

