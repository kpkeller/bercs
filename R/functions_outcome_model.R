# Functions to support fitting of the outcome models

# data -- dataframe of data
# Ht -- matrix/data frame of time splines to include
# return_sorted -- returns data frame that has been sorted, with added labels of
#                   'household_num
# log_transform -- log the concentrations?
##' create_standata_outcome()
##' @title Create List for Fitting Outcome Model via STAN
##' @description Generates a 'standata_outcome' object from a 'wide' data frame of outcome data
##' @param data data frame, or list of data frames, containg the data in 'wide' format. Expected to include (at least) "clust_id", "conc", "id", and "study"
##' @param xdf degrees of freedom to use in exposure spline.
##' @param Ht time spline matrix, or a list of such matrices.
##' @param ... arguments passed to \code{create_standata_outcome_singlestudy}
##' @seealso \code{\link{create_standata_exposure}}
##' @export
##' @importFrom Matrix bdiag
create_standata_outcome <- function(data, xdf=1, Ht=NULL,...) {
    if (length(data)>1 && class(data)=="list"){
        for (i in 1:length(data)){
            data[[i]]$clust_id <- paste0(data[[i]]$study, data[[i]]$clust_id)
            data[[i]]$id <- paste0(data[[i]]$study, data[[i]]$id)
        }
        data.new <- do.call(rbind, data)
        if (!is.null(Ht)){
            if (!all(sapply(Ht, inherits, "matrix"))) stop("'Ht' should be list of matrices.")
            for (i in 1:length(Ht)){
                class(Ht[[i]]) <- "matrix"
            }
            Ht.new <- Matrix::bdiag(Ht)
            Ht.new <- as.matrix(Ht.new)
        } else {
            Ht.new <- NULL
        }
    } else {
        data.new <- data
        Ht.new <- Ht
    }
    out <- create_standata_outcome_singlestudy(data=data.new, xdf=xdf, Ht=Ht.new,...)
    out
}

##' @rdname create_standata_outcome
##' @param covars character vector of variables to include as covariates, or a vector/matrix of covariate values. Should not include an intercept. If a character vector, this intercept is automatically removed (so default of "1" leads to no covariates).
##'
##' @param return_addition should the modified version of \code{data} be returned in addition to the \code{standata_outcome} object?
##' @param xfn name of function used to generate exposure splines
##' @param xfnargs named list of additional arguments for \code{xfn}
##' @param timefn name of function used to generate time splines
##' @param timedf degrees of freedom for time spline
##' @param timefnargs named list of additional arugments to \code{timefn}
##' @export
##' @importFrom splines2 iSpline
##' @importFrom stats model.matrix formula
create_standata_outcome_singlestudy <- function(data, xdf=1, Ht=NULL, covars="1",
                                                return_addition=FALSE,
                                                xfn="iSpline",xfnargs=list(),
                                                timefn="ns", timedf=0, timefnargs=list()){
    check_names(data, expected_names=c("clust_id", "conc", "id", "study"))
    check_names_to_overwrite(data, expected_names=c("subj_num", "cluster_num", "study_num"))

    out <- list()

    # Relabel ids into integer indicators
    # This facilitates linking back to original data and labels
    data$study_of_obs <- as.numeric(factor(data$study))
    data$cluster_of_obs <- as.numeric(factor(data$clust_id))
    data$subj_of_obs <- as.numeric(factor(data$id))

    out$S <- max(data$study_of_obs)
    out$K <- max(data$cluster_of_obs)
    out$n <- max(data$subj_of_obs)
    out$N <- nrow(data)
    out$study_of_obs <- data$study_of_obs
    out$subj_of_obs <- data$subj_of_obs
    out$cluster_of_subj <-  data$cluster_of_obs[!duplicated(data$subj_of_obs)][order(data$subj_of_obs[!duplicated(data$subj_of_obs)])]
    out$study_of_subj <- data$study_of_obs[!duplicated(data$subj_of_obs)][order(data$subj_of_obs[!duplicated(data$subj_of_obs)])]

    # Exposure values
    if (xdf==0){
        out$x <- 0
        Hx <- matrix(0, 0, 0)
    } else if (xdf==1){
        Hx <- scale(data$conc)
        out$x <- data$conc
    } else {
        out$x <- data$conc
        Hx <- do.call(create_spline_matrix,
                  c(list(x=data$conc,
                        df=xdf,
                        fn=xfn), xfnargs))
    }
    out$xdf <- xdf
    out$Hx <- as.matrix(Hx) # needs to be matrix class for STAN
    out$Hx_attributes <- attributes(Hx)

    # Covariates
    if (inherits(covars, "character")){
        Zformula <- paste0("~ ", paste0(covars, collapse="+"))
        Z <- stats::model.matrix(stats::formula(Zformula), data=data)
        Z <- Z[, -1, drop=FALSE] # Drop interecept
    } else if (inherits(covars, "vector")){
        Z <- as.matrix(Z)
    } else if (!inherits(covars, "matrix")){
        stop("'covars' must be character string, matrix, or vector.")
    }
    Z <- scale(Z)
    out$Z <- as.matrix(Z)
    if (ncol(out$Z)==0){
        out$Z <- matrix(0, 0, 0)
    }
    out$Zattr <- attributes(Z)
    out$p <- ncol(Z)

    # Outcome
    out$y <- data$case
    if (!is.null(data$at_risk)){
        out$nT <- data$at_risk
    } else {
        out$nT <- rep(1, nrow(data))
    }

    # Time
    if (is.null(Ht) && timedf>0){
        Ht <- do.call(create_spline_matrix,
                  c(list(x=data$date,
                        df=timedf,
                        fn=timefn), timefnargs))
    }
    if (!is.null(Ht)) {
        out <- add_Ht_standata(out, Ht)
    } else {
        out$Ht <- matrix(0, 0, 0)
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
                                 B=2000,
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
##' @param ... Passed to \code{\link{seq}} to control sequence of exposure values.
##'
##' @details This function creates a data frame containing the values of the exposure-response curve over a given range of values. The output is designed for easy plotting.
##' Currently, the fitted curve is plotted based upon the posterior means of the parameters, with the uncertainty bands based upon quantiles.
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
        bs_post <- cbind(bs_post, bs_post %*% table(standata$study_of_obs)/standata$N)
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
        if (!is.null(standata$Hx_attributes$`scaled:center`)) {
            exposure_seq_scaled <- (exposure_seq - standata$Hx_attributes$`scaled:center`)/standata$Hx_attributes$`scaled:scale`
        }
        else {
            exposure_seq_scaled <- exposure_seq
        }
    }
    else {
        Hxtemp <- standata$Hx
        attributes(Hxtemp) <- standata$Hx_attr
        predfn <- utils::getS3method(f = "predict",
                                     class(Hxtemp)[class(Hxtemp) !="matrix"][1])
        exposure_seq_scaled <- predfn(Hxtemp, newx = exposure_seq)
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
        if (nSout > nS) {
            bs_post <- cbind(bs_post, bs_post %*% table(standata$study_of_obs)/standata$N)
        }
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
##' @param incS Which studies to include. Default of NULL plots all studies.
##' @param ylab String providing y-axis label.
##' @param xlab String providing x-axis label.
##' @export
##' @import ggplot2
plot_fitted_ERC <- function (obj, incS=NULL, ylab = "Relative Risk", xlab = "Exposure")
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

    if(!is.null(incS)){
        fulldf <- subset(fulldf, study %in% incS)
    }
    ggplot(fulldf) + theme_bw() +
        geom_ribbon(aes(x = exposure,
                        ymin = exp(low),
                        ymax = exp(high),
                        group=factor(study)), fill = "grey80") +
        geom_line(aes(x = exposure,
                      y = exp(mean),
                      group=factor(study),
                      col=factor(study)),
                  lwd = 1.5) +
        geom_hline(yintercept = 1,
                   lty = 2) + xlab(xlab) + ylab(ylab)
}










