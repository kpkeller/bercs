# ########################
# Utility functions
#
# create_spline_matrix()
# create_spline()
# predict_spline()
# add_spline_time()
# add_spline_exposure()
# add_duplicated_flag()
#
# Unexported:
# check_names()
# check_names_to_overwrite()
# check_all_nonnegative_value()
# check_all_postive_value()
# ##############################



##' @title Create spline matrix
##' @description Wrapper function for creating a spline matrix
##' @param x Variable values for which spline is to be created.
##' @param xrange Range (min and max) of values to use for creating spline. Only used if \code{x} not provided.
##' @param by Interval of generated sequence.
##' @param df Degrees of freedom for spline. If set to 1, then a linear function assumed; otherwise \code{df} is passed to \code{fn}.
##' @param fn Character string specifying the spline-generating function. Defaults to \code{\link[splines]{ns}}.
##' @param ... additional arguments passed to \code{fn}.
##' @return For \code{create_spline}, a list containing:
##' \itemize{
##' \item{\code{x} -- the values provided in argument \code{x}}
##' \item{\code{spline} -- a matrix of spline values evaluated at \code{x}}
##' \item{\code{spline_mean} -- a vector of column-wise means of \code{spline}}
##' }
##' @export
##' @examples
##' x_ns <- create_spline(x=1:10, df=3, fn="ns")
##' x2_ns <- predict(x_ns$spline, newx=seq(0.5, 9.5, by=1))
##' x2_ns
##' @seealso \code{\link{add_spline_time}}, \code{\link{sim_update_time_splines}}
create_spline <- function(x, xrange=NULL, by=1, df=1, fn="ns",...){
    if(is.null(x)) {
        if (is.null(x)){
            stop("Either 'x' or 'xrange' required.")
        } else {
            x <- seq(xrange[1], to=xrange[2], by=by)
        }
    }
    spline <- create_spline_matrix(x, df=df, fn=fn, ...)
    out <- list(x=x,
                spline=spline,
                spline_mean=colMeans(spline))
    out
}

#' @rdname create_spline
##' @return For \code{create_spline_matrix}, the matrix of spline values only.
#' @export
create_spline_matrix <- function(x, df=1, fn="ns",...){
    if (!is.numeric(df) || df<1) stop("'df' must be postive integer.")
    if (df==1){
        Mx <- matrix(x)
    } else {
        fn <- get(fn)
        Mx <- fn(x=x,
                 df=df,
                 ...)
    }
    Mx
}

##' @rdname create_spline
##' @param spline Spline object from which to predict.
##' @param newx New values for which splines should be generated
##' @param center,scale Should spline values be centered and/or scaled? Default is \code{FALSE}. If a vector is provided, it is pased to \code{\link{scale}}.
##' @details Note that \code{predict_spline} relies on the existence of a \code{predict} generic for the class of spline.
##' @return For \code{predict_spline}, a matrix of spline values evaluated at \code{newx}.
##' @export
##' @importFrom utils getS3method
predict_spline <- function(spline, newx, center=FALSE, scale=FALSE){
    predfn <- utils::getS3method(f="predict", class(spline))
    H <- predfn(spline, newx=newx)

    if(any(as.logical(center)) || any(as.logical(scale))){
        H <- scale(H, center=center,scale=scale)
    }
    H
}


##' @title Add spline to a standata object
##' @description Updates a standata object with a matrix of time or exposure splines
##' @param standata the standata object to be updated
##' @param Mt,Mx the matrix of time or exposure splines to be added
##' @details \code{add_spline_time} updates the \code{Mt}, \code{timedf}, and \code{Mt_attributes} elements of \code{standata}. \code{add_spline_exposure} updates the \code{Mx}, \code{xdf}, and \code{Mx_attributes} elements of \code{standata}.
##' @return The \code{standata} object with updated elements.
##' @seealso \code{\link{create_standata_exposure}}, \code{\link{create_standata_outcome}}, \code{\link{expsim_update_parameter}}, \code{\link{outsim_update_parameter}}, \code{\link{create_spline}}, \code{\link{create_spline_matrix}}
#' @export
#' @examples
#' exp_data <- create_standata_exposure(group=rep(1, 10),
#'                                      conc=rnorm(10),
#'                                      unit_id=rep(0:1, 5),
#'                                      time=runif(10))
#' exp_data$Mt
#' Mt <- create_spline_matrix(exp_data$time, df=3, fn="ns")
#' exp_data <- add_spline_time(exp_data, Mt=Mt)
#' exp_data$Mt
add_spline_time <- function(standata, Mt){
    if(!inherits(Mt, c("matrix","numeric", "integer"))) stop("'Mt' should be a matrix or vector.")
    standata$Mt <- as.matrix(Mt)
    standata$timedf <- ncol(standata$Mt)
    standata$Mt_attributes <- attributes(Mt)
    standata
}

#' @rdname add_spline_time
#' @export
add_spline_exposure <- function(standata, Mx){
    if(!inherits(Mx, c("matrix","numeric", "integer"))) stop("'Mx' should be a matrix or vector.")
    standata$Mx <- as.matrix(Mx)
    standata$xdf <- ncol(standata$Mx)
    standata$Mx_attributes <- attributes(Mx)
    standata
}

#' @title Add Duplicated Observations Flag
#' @description Adds a column indicating duplicated rows
#' @param data Data frame to check.
#' @param dupvars Character vector providing the names of variables for which duplicated observations should be checked.
#' @details This is a wrapper around \code{\link{duplicated}}.
#' @return The data frame \code{data}, with an additional column \code{"dupobs"} that is a logical indicator of a duplicated row.
##' @export
add_duplicated_flag <- function(data, dupvars){
    check_names_to_overwrite(data, expected_names=c("dupobs"))
    data$dupobs <- duplicated(data[, dupvars])
    data
}


# Auxiliary functions to check names in data frame
# and either stop or warn if missing or will be modified, respectively
check_names <- function(df, expected_names){
    present <- expected_names %in% names(df)
    if(any(!present)){
        pp <- expected_names[!present]
        for (p in 1:length(pp)){
            stop(paste0("Column `", pp[p], "' required but not present."))
        }
    }
}

check_names_to_overwrite <- function(df, expected_names){
    present <- expected_names %in% names(df)
    if(any(present)){
        pp <- expected_names[present]
        for (p in 1:length(pp)){
            warning(paste0("Column `", pp[p], "' will be modified"))
        }
    }
}

check_all_nonnegative_value <- function(x){
    all(is.numeric(x)) && all(x>=0) && any(x>0)
}
check_all_postive_value <- function(x){
    all(is.numeric(x)) && all(x>0)
}

