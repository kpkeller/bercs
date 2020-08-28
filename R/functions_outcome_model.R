# Functions to support fitting of the outcome models

##' @title Create List for Fitting Outcome Model via STAN
##' @description Generates a 'standata_outcome' object from a 'long' data frame of outcome data
##' @param ... arguments passed to \code{create_standata_outcome_singlestudy}. See Details.
##' @param datalist List of data frames, containing the data in 'long' format.
##' @param Mtlist time spline matrix, or a list of such matrices.
##' @param covarslist list of covariate variable names for each study, or a list of covariate matrices.
##' @details An important special case of \code{create_standata_outcome} is when there is only a single study. In this case, the values are passed to \code{create_standata_outcome_singlestudy} without pre-processing.
##'
##' If using a separate temporal spline for each study, first compute the time splines for each study individually using \code{\link{create_spline}}. Combine these into a list and pass that to \code{Mtlist}.
##'
##' The matrices of adjustment variables are created using the names provided to \code{covarslist}. This is done using \code{\link[stats]{formula}} and \code{\link[stats]{model.matrix}}, so transformations and interactions can be used in the typical manner.
##' @seealso \code{\link{create_standata_exposure}}, \code{\link{add_spline_exposure}}, \code{\link{sample_outcome_model}}
##' @export
##' @importFrom Matrix bdiag
create_standata_outcome <- function(...,
                                    datalist=NULL,
                                    Mtlist=NULL,
                                    covarslist=NULL) {
    if (!is.null(datalist)){
        if (length(datalist)>1 && class(datalist)=="list"){
            for (i in 1:length(datalist)){
              if (!all(is.null(datalist[[i]]$clust_id))){
                datalist[[i]]$clust_id <- paste0(datalist[[i]]$study, datalist[[i]]$clust_id)
              }
                datalist[[i]]$unit_id <- paste0(datalist[[i]]$study, datalist[[i]]$unit_id)
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
            if (!is.null(covarslist)){
              covarslist.new <- vector("list", length(covarslist))
              for (i in 1:length(covarslist)){
                if (inherits(covarslist[[i]], "character")){
                  Zformula <- paste0("~ ", paste0(covarslist[[i]], collapse="+"))
                  Z <- stats::model.matrix(stats::formula(Zformula), data=datalist[[i]])
                  Z <- Z[, -1, drop=FALSE] # Drop intercept
                } else if (inherits(covarslist[[i]], "vector")){
                  Z <- as.matrix(covarslist[[i]])
                } else if (inherits(covarslist[[i]], "matrix")){
                  Z <- covarslist[[i]]
                } else {
                  stop("'covarslist' must be a list containing a character string, matrix, or vector.")
                }
                covarslist.new[[i]] <- Z
              }
              covars.new <- Matrix::bdiag(covarslist.new)
              covars.new <- as.matrix(covars.new)
              colnames(covars.new) <- unlist(sapply(covarslist.new, colnames))
            } else {
              covars.new <- c("1")
            }
        } else {
            data.new <- datalist
            Mt.new <- Mtlist
        }
        out <- create_standata_outcome_singlestudy(data=data.new,
                                                   Mt=Mt.new,
                                                   covars=covars.new,
                                                   ...)
    } else {
        out <- create_standata_outcome_singlestudy(...)
    }
    out
}

##' @rdname create_standata_outcome
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
        Mx <- scale(newdata$conc)
        out$x <- newdata$conc
    } else {
        out$x <- newdata$conc
        Mx <- do.call(create_spline_matrix,
                  c(list(x=newdata$conc,
                        df=xdf,
                        fn=xfn), xfnargs))
    }
    out <- add_spline_exposure(out, Mx=Mx)

    # Covariates
    if (inherits(covars, "character")){
        Zformula <- paste0("~ ", paste0(covars, collapse="+"))
        Z <- stats::model.matrix(stats::formula(Zformula), data=data)
        Z <- Z[, -1, drop=FALSE] # Drop intercept
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
##' @examples
##' data(casedataA)
##' outcome_dataA <- create_standata_outcome(data=casedataA)
##' outcome_dataA <- add_priors(outcome_dataA,
##'                             sigmaI=c(0, 0.1))
##' outcome_mod_fit1 <- sample_outcome_model(outcome_dataA)
##' print(outcome_mod_fit1, pars=c("reI_raw", "reI","mui"), include=FALSE)
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
##' @param exprange range of exposure values over which to compute the curve
##' @param expsequence sequence of exposure values at which the curve should be evaluated. For plotting, it is preferable to use \code{exprange}, but for specific known exposure values use \code{expsequence}.
##' @param ciband width of credible interval band
##' @param inclInterceptUncertainty Should intercept uncertainty be included uncertainty estimates? See details for more information.
##' @param inclIntercept Should the intercept term be included in the curve?
##' @param intercept_prop Proportions used in calculating the "average" intercept. Defaults to equal proportions for each study (\code{"equal"}) and can be set to be proportional to the number of observations in each study (\code{"obs"}). A vector of proportions can also be given.
##' @param study Return the curve for a specific study, or the
##' @param beta_post vector or matrix of posterior samples for 'beta' parameter. Only needed if \code{stanfit} not provided.
##' @param bs_post vector or matrix of posterior samples for 'bS' parameter. Only needed if \code{stanfit} not provided.
##' @param nS Number of studies. Only needed if \code{standata} not provided.
##' @param Mx Spline matrix for exposure. Only needed if \code{standata} not provided.
##' @param Mx_attributes Attributes for the spline matrix. Only needed if \code{standata} not provided.
##' @param xdf Degrees of freedom in exposure splines
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
compute_ERC <- function (standata,
                            stanfit,
                            exprange = c(0, 100),
                            expsequence=NULL,
                            ref_exposure=NULL,
                            ciband = 0.95,
                            inclInterceptUncertainty=TRUE,
                            inclIntercept=FALSE,
                            intercept_prop=c("equal", "obs"),
                            study=NULL,
                            beta_post,
                            bs_post,
                            nS=standata$S,
                            Mx=standata$Mx,
                            Mx_attributes=standata$Mx_attributes,
                            xdf=standata$xdf,
                            ...)
{
    if (ciband < 0 || ciband > 1)
        stop("'ciband' must be between 0 and 1.")
    if (missing(beta_post)) {
        beta_post <- rstan::extract(stanfit, pars = "beta")$beta
    }
    # If multiple studies, will add column with "average" intercept
    nSout <- ifelse(nS>1 & length(dim(beta_post))==2, nS+1, nS)
    if (is.null(expsequence)){
      expsequence <- seq(exprange[1], exprange[2],...)
    } else {
      if (!is.null(ref_exposure)){
        if (!ref_exposure %in% expsequence) {
          expsequence <- c(ref_exposure, expsequence)
        }
      }
    }
    if (missing(bs_post)) {
      bs_post <- extract(stanfit, pars = "bS")$bS
    }

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

    beta_post <- reshape_beta_post(beta_post=beta_post,
                      nSout=nSout)

    # if (length(dim(beta_post))==2) {
    #     beta_post2 <- array(dim=c(nrow(beta_post),
    #                               nSout,
    #                               ncol(beta_post)))
    #     for (j in 1:dim(beta_post2)[2]){
    #         beta_post2[, j, ] <-beta_post
    #     }
    #     beta_post <- beta_post2
    # }



        #
    # if (xdf == 1) {
    #     if (!is.null(Mx_attributes$`scaled:center`)) {
    #         exposure_seq_scaled <- (expsequence - Mx_attributes$`scaled:center`)/Mx_attributes$`scaled:scale`
    #     }
    #     else {
    #         exposure_seq_scaled <- expsequence
    #     }
    # }
    # else {
    #     Mxtemp <- Mx
    #     attributes(Mxtemp) <- Mx_attributes
    #     predfn <- utils::getS3method(f = "predict",
    #                                  class(Mxtemp)[class(Mxtemp) !="matrix"][1])
    #     exposure_seq_scaled <- predfn(Mxtemp, newx = expsequence)
    # }
    # fitted_seq <- array(dim=c(length(expsequence),
    #                           dim(beta_post)[1], # B from stan fit
    #                           nSout))
    # for (i in 1:nSout){
    #     fitted_seq[, , i] <- exposure_seq_scaled %*% t(beta_post[, i, ])
    # }

    fitted_seq <- compute_fitted_sequence(beta_post = beta_post,
                                          expsequence=expsequence,
                                          nSout=nSout,
                                          xdf=xdf,
                                          Mx=Mx,
                                          Mx_attributes=Mx_attributes)

    if (inclIntercept & !inclInterceptUncertainty) {
        inclInterceptUncertainty <- TRUE
        warning("inclIntercept set to TRUE. This requires inclInterceptUncertainty to be TRUE.")
    }
    if (inclInterceptUncertainty) {
        if (!inclIntercept){
            bs_post <- sweep(bs_post, 2, colMeans(bs_post), FUN = "-")
        }
        bs_post <- rep(1, length(expsequence)) %o% bs_post

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
                exposure = expsequence)
    dflist <- vector("list", nSout)
    for (j in 1:nSout){
      dflist[[j]] <- data.frame(mean=obj$mean[, j],
                                low=obj$low[, j],
                                high=obj$high[, j],
                                exposure=obj$exposure,
                                study=j)
    }
    if (length(dflist)==1){
      dflist <- dflist[[1]]
    }
    if (!is.null(ref_exposure)){
      dflist <- center_ERC(dflist,
                        ref_exposure=ref_exposure)
    }
    dflist
}

##' @rdname compute_ERC
##' @param obj Data frame containing \code{exposure}, \code{mean}, \code{low}, and \code{high}. Typically generated from \code{\link{compute_ERC}}.
##' @param incS If model has a single curve but with different selections of intercept uncertainty, this indicates which choice of uncertainty to use. Defaults to the last column of the curve in \code{obj}, which is typically the averaged intercept uncertainty. If multiple curves are fit, this selects which curve(s) is plotted.
##' @param expERC Should the fitted curve be exponentiated (TRUE) or not (FALSE).
##' @param ylab String providing y-axis label.
##' @param xlab String providing x-axis label.
##' @param ribbon Should the uncertainty be represented as a filled ribbon (TRUE) or lines without fill (FALSE).
##' @export
##' @import ggplot2
plot_ERC <- function (obj,
                             incS=NULL,
                             expERC=TRUE,
                             ylab = "Relative Risk",
                             xlab = "Exposure",
                             ribbon=FALSE)
{
    if (is.list(obj)){
      nS <- length(obj)
    } else {
      nS <- 1
    }
    fulldf <- do.call(rbind, obj)
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
        g <- g + geom_ribbon(aes(x = .data$exposure,
                                 ymin = .data$low,
                                 ymax = .data$high,
                                 group=.data$study),
                             fill = "grey80")
    } else {
        g <- g+   geom_line(aes(x=.data$exposure,
                                y=.data$low,
                                group = .data$study,
                                col=.data$study),
                            data=fulldf) +
            geom_line(aes(x=.data$exposure,
                          y=.data$high,
                          group = .data$study,
                          col=.data$study),
                      data=fulldf)
    }
      g   + geom_line(aes(x = .data$exposure,
                      y = .data$mean,
                      group=.data$study,
                      col=.data$study),
                  lwd = 1.5) +
        geom_hline(yintercept = 1,
                   lty = 2) + xlab(xlab) + ylab(ylab)

}


##' @rdname compute_ERC
##' @param ref_exposure Exposure concentration at which the spline should be shifted to have value 0 on log scale (value 1 on exponentiated scale). The largest exposure value less than or equal to \code{ref_exposure} is used as the new reference concentration.
##' @export
center_ERC <- function(obj, ref_exposure=min(obj$exposure)){
    obj_new <- obj
    if (is.list(obj_new)){
      for (j in 1:length(obj_new)){
        ind_cut <- max(which(obj_new[[j]]$exposure <= ref_exposure))
        obj_new[[j]]$low <- obj_new[[j]]$low - obj_new[[j]]$mean[ind_cut]
        obj_new[[j]]$high <- obj_new[[j]]$high - obj_new[[j]]$mean[ind_cut]
        obj_new[[j]]$mean <- obj_new[[j]]$mean - obj_new[[j]]$mean[ind_cut] # do last
      }} else {
        ind_cut <- max(which(obj_new$exposure <= ref_exposure))
        obj_new$low <- obj_new$low - obj_new$mean[ind_cut]
        obj_new$high <- obj_new$high - obj_new$mean[ind_cut]
        obj_new$mean <- obj_new$mean - obj_new$mean[ind_cut] # do last
      }
    obj_new
}



##' @rdname compute_ERC
##' @details To calculate odds ratios for specific combinations of exposure values, use the `compute_OR` function, which will correctly calculate the credible intervals for the relative difference.
##' @export
compute_OR <- function (standata,
                        stanfit,
                        exprange = c(0, 100),
                        expsequence=NULL,
                        ref_exposure=0,
                        ciband = 0.95,
                        inclInterceptUncertainty=TRUE,
                        inclIntercept=FALSE,
                        intercept_prop=c("equal", "obs"),
                        study=NULL,
                        beta_post,
                        bs_post,
                        nS=standata$S,
                        Mx=standata$Mx,
                        Mx_attributes=standata$Mx_attributes,
                        xdf=standata$xdf,
                        ...)
{
  if (ciband < 0 || ciband > 1)
    stop("'ciband' must be between 0 and 1.")
  if (missing(beta_post)) {
    beta_post <- rstan::extract(stanfit, pars = "beta")$beta
  }

  if (is.null(expsequence)){
    expsequence <- seq(exprange[1], exprange[2],...)
  } else {
    if (!is.null(ref_exposure)){
      if (!ref_exposure %in% expsequence) {
        expsequence <- c(ref_exposure, expsequence)
      }
    }
  }

  # Are there multiple curves? Usually no.
  multiple_curves <- ifelse(length(dim(beta_post))>2, TRUE, FALSE)

  beta_post <- reshape_beta_post(beta_post,
                    nSout=nS)
  fitted_seq <- compute_fitted_sequence(beta_post = beta_post,
                                        expsequence=expsequence,
                                        nSout=nS,
                                        xdf=xdf,
                                        Mx=Mx,
                                        Mx_attributes=Mx_attributes)

  fitted_seq_ref <- fitted_seq[which(expsequence==ref_exposure),,]
  fitted_seq <- sweep(fitted_seq, 2:3, fitted_seq_ref,FUN="-")
  fitted_seq_mean <- apply(fitted_seq, c(1, 3), mean)
  fitted_seq_low <- apply(fitted_seq, c(1,3), stats::quantile,
                          probs = (1 - ciband)/2)
  fitted_seq_high <- apply(fitted_seq, c(1,3), stats::quantile,
                           probs = 0.5 + ciband/2)

  obj <- list(exposure = expsequence,
              logOR_mean = fitted_seq_mean,
              logOR_low = fitted_seq_low,
              logOR_high = fitted_seq_high,
              OR_mean = exp(fitted_seq_mean),
              OR_low = exp(fitted_seq_low),
              OR_high = exp(fitted_seq_high))
  dflist <- vector("list", nS)
  for (j in 1:nS){
    dflist[[j]] <- data.frame(exposure=obj$exposure,
                              logOR_mean=obj$logOR_mean[, j],
                              logOR_low=obj$logOR_low[, j],
                              logOR_high=obj$logOR_high[, j],
                              OR_mean=obj$OR_mean[, j],
                              OR_low=obj$OR_low[, j],
                              OR_high=obj$OR_high[, j],
                              study=j)
  }
  if (!multiple_curves){
    dflist <- dflist[[1]]
  }
  dflist
}


reshape_beta_post <- function(beta_post, nSout){
  if (length(dim(beta_post))==2) {
    beta_post2 <- array(dim=c(nrow(beta_post),
                              nSout,
                              ncol(beta_post)))
    for (j in 1:dim(beta_post2)[2]){
      beta_post2[, j, ] <-beta_post
    }
    beta_post <- beta_post2
  }
  beta_post
}

# Compute the fitted sequence from an exposure sequence
# and set of beta values
## @param beta_post Matrix of posterior sample of beta's
## @param expsequence Sequence of exposure values
## @param nSout Number of studies + 1 (unless only one study) in model output
## @param xdf Degrees of freedom in exposure spline
## @param Mx Matrix of exposure spline values
## @param Mx_attributes Attributes of Mx
compute_fitted_sequence <- function(beta_post,
                                    expsequence,
                                    nSout,
                                    xdf,
                                    Mx,
                                    Mx_attributes){
if (xdf == 1) {
  if (!is.null(Mx_attributes$`scaled:center`)) {
    exposure_seq_scaled <- (expsequence - Mx_attributes$`scaled:center`)/Mx_attributes$`scaled:scale`
  }
  else {
    exposure_seq_scaled <- expsequence
  }
}
else {
  Mxtemp <- Mx
  attributes(Mxtemp) <- Mx_attributes
  predfn <- utils::getS3method(f = "predict",
                               class(Mxtemp)[class(Mxtemp) !="matrix"][1])
  exposure_seq_scaled <- predfn(Mxtemp, newx = expsequence)
}
fitted_seq <- array(dim=c(length(expsequence),
                          dim(beta_post)[1], # B from stan fit
                          nSout))
for (i in 1:nSout){
  fitted_seq[, , i] <- exposure_seq_scaled %*% t(beta_post[, i, ])
}
  fitted_seq
}
