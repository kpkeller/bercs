# Functions to support fitting of the outcome models


# return_sorted -- returns data frame that has been sorted, with added labels of
#                   'household_num
# log_transform -- log the concentrations?
##' create_standata_outcome()
##' @title Create List for Fitting Outcome Model via STAN
##' @description Generates a 'standata_outcome' object from a 'long' data frame of outcome data
##' @param datalist List of data frames, containing the data in 'long' format.
##' @param Mtlist time spline matrix, or a list of such matrices.
##' @param ... arguments passed to \code{create_standata_outcome_singlestudy}. See Details.
##' @details An important special case of \code{create_standata_outcome} is when there is only a single study. In this case, all arguments should be named so that they are correctly passed to \code{create_standata_outcome_singlestudy} without pre-processing.
##' @seealso \code{\link{create_standata_exposure}}
##' @export
##' @importFrom Matrix bdiag
create_standata_outcome <- function(datalist=NULL,
                                    Mtlist=NULL,
                                    ...) {
    if (!is.null(datalist)){
        if (length(datalist)>1 && class(datalist)=="list"){
            for (i in 1:length(datalist)){
                datalist[[i]]$clust_id <- paste0(datalist[[i]]$study, datalist[[i]]$clust_id)
                datalist[[i]]$id <- paste0(datalist[[i]]$study, datalist[[i]]$id)
            }
            data.new <- do.call(rbind, datalist)
            if (!is.null(Mtlist)){
                if (!all(sapply(Mtlist, inherits, "matrix"))) stop("'Mtlist' should be list of matrices.")
                # for (i in 1:length(Mtlist)){
                #     class(Mtlist[[i]]) <- "matrix"
                # }
                Mt.new <- Matrix::bdiag(Mtlist)
                Mt.new <- as.matrix(Mt.new)
            } else {
                Mt.new <- NULL
            }
        } else {
            data.new <- datalist
            Mt.new <- Mtlist
        }
        out <- create_standata_outcome_singlestudy(data=data.new, Mt=Mt.new,...)
    } else {
        out <- create_standata_outcome_singlestudy(...)
    }
    out
}

##' @title create_standata_outcome
##'
##' @param return_addition should the modified version of \code{data} be returned in addition to the \code{standata_outcome} object?
##' @param data Optional data frame containing data in a long format.
##' @param unit_id Vector identifying the distinct unit (e.g. person) for each observation.
##' @param conc Vector of exposure concentrations for each observation.
##' @param study Vector identifying the study of each observation.
##' @param clust_id Optional vector identifying cluster membership for each observation.
##' @param case Count of cases for the observations. In most cases, this is a 0/1 indicator.
##' @param at_risk Time at risk for the observation
##' @param covars character vector of variables to include as covariates, or a vector/matrix of covariate values. Should not include an intercept. If a character vector, this intercept is automatically removed (so default of "1" leads to no covariates).
##' @param scale_covars should the covariate matrix be centered and scaled?
##' @param xdf Degrees of freedom for exposure splines
##' @param Mt Optional matrix of time splines
##' @param xfn name of function used to generate exposure splines
##' @param xfnargs named list of additional arguments for \code{xfn}
##' @param timefn name of function used to generate time splines
##' @param timedf degrees of freedom for time spline
##' @param timefnargs named list of additional arugments to \code{timefn}
##' @export
##' @importFrom splines2 iSpline
##' @importFrom stats model.matrix formula
create_standata_outcome_singlestudy <- function(data=NULL,
                                                unit_id=data$unit_id,
                                                conc=data$conc,
                                                study=data$study,
                                                clust_id=data$clust_id,
                                                case=data$case,
                                                at_risk=data$at_risk,
                                                covars="1",
                                                scale_covars=TRUE,
                                                xdf=1,
                                                xfn="iSpline",
                                                xfnargs=list(),
                                                Mt=NULL,
                                                timefn="ns",
                                                timedf=0,
                                                timefnargs=list(),
                                                return_addition=FALSE){

    N <- length(conc)
    if(!all(length(unit_id)==N,
            length(case)==N)) stop("The lengths of 'unit_id', 'case', and 'conc' must be equal")
    if (xdf < 0) stop("'xdf' must be a non-negative integer.")
    # Create data frame
    # Relabel ids into integer indicators
    # This facilitates linking back to original data and labels
    newdata <- data.frame(conc=conc,
                          unit_of_obs=as.numeric(factor(unit_id)),
                          case=case)
    if (is.null(at_risk)){
        newdata$at_risk <- rep(1, N)
    } else {
        newdata$at_risk <- at_risk
    }
    if (!is.null(clust_id) && any(clust_id!=0)){
        if (length(clust_id)!=N) stop("The lengths of 'clust_id' and 'conc' must be equal")
        nocluster <- FALSE
        newdata$cluster_of_obs <- as.numeric(factor(clust_id))
    } else {
        nocluster <- TRUE
        newdata$cluster_of_obs <- rep(0, N)
    }
    if (!is.null(study)){
        if (length(study)!=N) stop("The lengths of 'study' and 'conc' must be equal")
        newdata$study_of_obs <- as.numeric(factor(study))
    } else {
        newdata$study_of_obs <- rep(1, N)
    }

    out <- list()
    out$S <- max(newdata$study_of_obs)
    out$K <- max(newdata$cluster_of_obs)
    out$n <- max(newdata$unit_of_obs)
    out$N <- N
    out$study_of_obs <- newdata$study_of_obs
    out$unit_of_obs <- newdata$unit_of_obs
    out$cluster_of_unit <-  newdata$cluster_of_obs[!duplicated(newdata$unit_of_obs)][order(newdata$unit_of_obs[!duplicated(newdata$unit_of_obs)])]
    out$study_of_unit <- newdata$study_of_obs[!duplicated(newdata$unit_of_obs)][order(newdata$unit_of_obs[!duplicated(newdata$unit_of_obs)])]

    # Exposure values
    if (xdf==0){
        out$x <- 0
        Mx <- matrix(0, 0, 0)
    } else if (xdf==1){
        Mx <- scale(data$conc)
        out$x <- data$conc
    } else {
        out$x <- data$conc
        Mx <- do.call(create_spline_matrix,
                  c(list(x=data$conc,
                        df=xdf,
                        fn=xfn), xfnargs))
    }
    out <- add_spline_exposure(out, Mx=Mx)

    # Covariates
    if (inherits(covars, "character")){
        Zformula <- paste0("~ ", paste0(covars, collapse="+"))
        Z <- stats::model.matrix(stats::formula(Zformula), data=data)
        Z <- Z[, -1, drop=FALSE] # Drop interecept
    } else if (inherits(covars, "vector")){
        Z <- as.matrix(covars)
    } else if (inherits(covars, "matrix")){
        Z <- covars
    } else {
        stop("'covars' must be character string, matrix, or vector.")
    }
    if (scale_covars){
      Z <- scale(Z)
    }
    out$Z <- as.matrix(Z)
    if (ncol(out$Z)==0){
        out$Z <- matrix(0, 0, 0)
    }
    out$Zattr <- attributes(Z)
    out$p <- ncol(Z)

    # Outcome
    out$y <- newdata$case
    out$nT <- newdata$at_risk

    # Time
    if (is.null(Mt) && timedf>0){
        Mt <- do.call(create_spline_matrix,
                  c(list(x=data$date,
                        df=timedf,
                        fn=timefn), timefnargs))
    }
    if (!is.null(Mt)) {
        out <- add_spline_time(out, Mt)
    } else {
        out$Mt <- matrix(0, 0, 0)
        out$timedf <- 0
    }

    class(out) <- "standata_outcome"
    if (return_addition){
        return(list(standata=out,
                    df=data))
    } else {
        return(out)
    }
}



##' @title Sample Outcome Model
##' @description Samples from the posterior distribution of outcome model parameters, using STAN
##' @param standata An object of class `standata_outcome`, typically created from \code{\link{create_standata_outcome}}.
#' @inheritParams sample_exposure_model
##' @param multiple_exposure_curves Logical indicating whether exposure response curves should be estimated separately by study.
##' @param restrictBeta Logical indicating that model should be fit that forces a positive value for the exposure spline coefficients.
##' @author Joshua Keller
##' @seealso \link{create_standata_outcome}, \link{outsim_sample_observations}
##' @export
sample_outcome_model <- function(standata,
                                 B=1000,
                                 warmup=B,
                                 chains=4,
                                 control=list(adapt_delta=0.9,
                                              max_treedepth=12),
                                 multiple_exposure_curves=FALSE,
                                 restrictBeta=FALSE,
                                 ...){

    if (!inherits(standata, "standata_outcome")) stop("`standata` must be of class 'standata_outcome'.")

    if(restrictBeta){
        if (!is.null(standata$beta_lower_lim) && !is.na(standata$beta_lower_lim)){
            warning("Because 'restrictBeta=TRUE', replacing existing value of 'standata$beta_lower_lim' with 0.")
        }
        standata$beta_lower_lim <- 0
    }
    if (is.numeric(standata$beta_lower_lim) && standata$xdf==0){
        stop("Cannot restrict beta when xdf==0. Set xdf to be >0 or beta_lower_lim to be NULL")
    }
    if (multiple_exposure_curves && standata$xdf==0){
        stop("Cannot have multiple exposure curves when xdf==0. Set xdf to be >0 or multiple_exposure_curves to be FALSE")
    }

    model_name <- "outcome_model"
    if (multiple_exposure_curves) model_name <- paste0(model_name, "_multipleExpCurves")
    if (is.numeric(standata$beta_lower_lim))  model_name <- paste0(model_name, "_restrictBeta")

    out_stanfit <- sampling(stanmodels[[model_name]],
                            data = standata,
                            iter = B + warmup,
                            warmup=warmup,
                            chains = chains,
                            control=control,
                            ...)

    out_stanfit
}



# Functions for post-processing outcome model results

##' @title Extract fitted Exposure-Response Curve (ERC)
##' @description Computes the exposure-response curve from a fitted outcome model
##' @param standata \code{standata_outcome} object used to fit the model
##' @param stanfit fitted outcome model from \code{\link{sample_outcome_model}}. Posterior samples of 'beta' are extracted from this object.
##' @param beta_post vector or matrix of posterior samples for 'beta' parameter. Only needed if \code{stanfit} not provided.
##' @param exprange range of exposure values over which to compute the curve
##' @param ciband width of credible interval band
##' @param inclInterceptUncertainty Should intercept uncertainty be included uncertainty estimates? See details for more information.
##' @param inclIntercept Should the intercept term be included in the curve?
##' @param intercept_prop Proportions used in calculating the "average" intercept. Defaults to equal proportions for each study (\code{"equal"}) and can be set to be proportional to the number of observations in each study (\code{"obs"}). A vector of proportions can also be given.
##' @param ... Passed to \code{\link{seq}} to control sequence of exposure values.
##'
##' @details This function creates a data frame containing the values of the exposure-response curve over a given range of values. The output is designed for easy plotting.
##' Currently, the fitted curve is plotted based upon the posterior means of the parameters, with the uncertainty bands based upon quantiles.
##'
##' Uncertainty from the intercepts is included in the confidence bands by default, since this corresponds to a common interpretation of such intervals. This requires picking a single value for the intercept. For models fit to multiple studies, a separate set of uncertainty will be created for each study. Additionally, a set of results corresponding to "average" intercept is created. The \code{intercept_prop} argument controls the relative contribution of the intercepts from each model.
##'
##' @seealso \code{\link{sample_outcome_model}}
##' @export
#' @importFrom stats quantile
#' @importFrom utils getS3method
#' @importFrom rstan extract
get_fitted_ERC <- function (standata,
                            stanfit,
                            beta_post,
                            exprange = c(0, 100),
                            ciband = 0.95,
                            inclInterceptUncertainty=TRUE,
                            inclIntercept=FALSE,
                            intercept_prop=c("equal", "obs"),
                            ...)
{
    nS <- standata$S
    if (ciband < 0 || ciband > 1)
        stop("'ciband' must be between 0 and 1.")
    if (missing(beta_post)) {
        beta_post <- extract(stanfit, pars = "beta")$beta
    }
    # If multiple studies, will add column with "average" intercept
    nSout <- ifelse(nS>1 & length(dim(beta_post))==2, nS+1, nS)
    exposure_seq <- seq(exprange[1], exprange[2],...)
    bs_post <- extract(stanfit, pars = "bS")$bS

    if (nSout > nS){
        # Add column that has "average" intercept, if there are multiple studies
        if (is.null(intercept_prop)){
            intercept_prop <- table(standata$study_of_obs)/standata$N
        }
        if (intercept_prop=="equal"){
            intercept_prop <- rep(1/nS, length=nS)
        } else if (intercept_prop=="obs"){
            intercept_prop <- table(standata$study_of_obs)/standata$N
        } else if (!is.numeric(intercept_prop)){
            stop("'intercept_prop' must be 'equal', 'open' or a numeric vector")
        }
        if (sum(intercept_prop)!=1) stop("The values in 'intercept_prop' should equal 1.")
        if (any(intercept_prop < 0)) stop("The values in 'intercept_prop' should be non-negative.")


        bs_post <- cbind(bs_post, bs_post %*% intercept_prop)
    }
    if (length(dim(beta_post))==2) {
        beta_post2 <- array(dim=c(nrow(beta_post),
                                  nSout,
                                  ncol(beta_post)))
        for (j in 1:dim(beta_post2)[2]){
            beta_post2[, j, ] <-beta_post
        }
        beta_post <- beta_post2
    }

    if (standata$xdf == 1) {
        if (!is.null(standata$Mx_attributes$`scaled:center`)) {
            exposure_seq_scaled <- (exposure_seq - standata$Mx_attributes$`scaled:center`)/standata$Mx_attributes$`scaled:scale`
        }
        else {
            exposure_seq_scaled <- exposure_seq
        }
    }
    else {
        Mxtemp <- standata$Mx
        attributes(Mxtemp) <- standata$Mx_attr
        predfn <- utils::getS3method(f = "predict",
                                     class(Mxtemp)[class(Mxtemp) !="matrix"][1])
        exposure_seq_scaled <- predfn(Mxtemp, newx = exposure_seq)
    }
    fitted_seq <- array(dim=c(length(exposure_seq),
                              dim(beta_post)[1], # B from stan fit
                              nSout))
    for (i in 1:nSout){
        fitted_seq[, , i] <- exposure_seq_scaled %*% t(beta_post[, i, ])
    }

    if (inclIntercept & !inclInterceptUncertainty) {
        inclInterceptUncertainty <- TRUE
        warning("inclIntercept set to TRUE. This requires inclInterceptUncertainty to be TRUE.")
    }
    if (inclInterceptUncertainty) {
        if (!inclIntercept){
            bs_post <- sweep(bs_post, 2, colMeans(bs_post), FUN = "-")
        }
        bs_post <- rep(1, length(exposure_seq)) %o% bs_post

        fitted_seq <- fitted_seq + bs_post
    }
    fitted_seq_mean <- apply(fitted_seq, c(1, 3), mean)
    fitted_seq_low <- apply(fitted_seq, c(1,3), stats::quantile,
                                  probs = (1 - ciband)/2)
    fitted_seq_high <- apply(fitted_seq, c(1,3), stats::quantile,
                                   probs = 0.5 + ciband/2)

    obj <- list(mean = fitted_seq_mean,
                low = fitted_seq_low,
                high = fitted_seq_high,
                exposure = exposure_seq)
    obj
}

##' @rdname get_fitted_ERC
##' @param obj Data frame containing \code{exposure}, \code{mean}, \code{low}, and \code{high}. Typically generated from \code{\link{get_fitted_ERC}}.
##' @param incS If model has a single curve but with different selections of intercept uncertainty, this indicates which choice of uncertainty to use. Defaults to the last column of the curve in \code{obj}, which is typically the averaged intercept uncertainty. If multiple curves are fit, this selects which curve(s) is plotted.
##' @param expERC Should the fitted curve be exponentiated (TRUE) or not (FALSE).
##' @param ylab String providing y-axis label.
##' @param xlab String providing x-axis label.
##' @param ribbon Should the uncertainty be represented as a filled ribbon (TRUE) or lines without fill (FALSE).
##' @export
##' @import ggplot2
plot_fitted_ERC <- function (obj, incS=NULL, expERC=TRUE, ylab = "Relative Risk", xlab = "Exposure", ribbon=FALSE)
{
    nS <- ncol(obj$mean)
    dflist <- vector("list", nS)
    for (j in 1:nS){
        dflist[[j]] <- data.frame(mean=obj$mean[, j],
                                  low=obj$low[, j],
                                  high=obj$high[, j],
                                  exposure=obj$exposure,
                                  study=j)
    }
    fulldf <- do.call(rbind, dflist)
    if (expERC) {
        fulldf$mean <- exp(fulldf$mean)
        fulldf$low <- exp(fulldf$low)
        fulldf$high <- exp(fulldf$high)
    }

    if(is.null(incS)){
        incS <- nS
    }
    fulldf <- subset(fulldf, study %in% incS)
    fulldf$study <- factor(fulldf$study)
    g <- ggplot(fulldf) + theme_bw()
    if (ribbon){
        g <- g + geom_ribbon(aes(x = exposure,
                                 ymin = low,
                                 ymax = high,
                                 group=study),
                             fill = "grey80")
    } else {
        g <- g+   geom_line(aes(x=exposure,
                                y=low,
                                group = study,
                                col=study),
                            data=fulldf) +
            geom_line(aes(x=exposure,
                          y=high,
                          group = study,
                          col=study),
                      data=fulldf)
    }
      g   + geom_line(aes(x = exposure,
                      y = mean,
                      group=study,
                      col=study),
                  lwd = 1.5) +
        geom_hline(yintercept = 1,
                   lty = 2) + xlab(xlab) + ylab(ylab)

}


##' @rdname get_fitted_ERC
##' @param ref_exposure Exposure concentration at which the spline should be shifted to have value 0 on log scale (value 1 on exponentiated scale). The largest exposure value less than \code{ref_exposure} is used as the new reference concentration.
##' @export
center_ERC <- function(obj, ref_exposure=min(obj$exposure)){
    obj_new <- obj
    ind_cut <- max(which(obj$exposure < ref_exposure))
    for (j in 1:ncol(obj_new$mean)){
        obj_new$mean[, j] <- obj$mean[, j] - obj$mean[ind_cut, j]
        obj_new$low[, j] <- obj$low[, j] - obj$mean[ind_cut, j]
        obj_new$high[, j] <- obj$high[, j] - obj$mean[ind_cut, j]
    }
    obj_new
}






