##' @title Create spline matrix
##' @description Wrapper function for creating a spline matrix
##' @param x Variable values for which spline is to be created.
##' @param df Degrees of freedom for spline. If set to 1, then a linear function assumed; otherwise \code{df} is passed to \code{fn}.
##' @param fn Character string specifying the spline-generating function. Defaults to \code{\link[splines]{ns}}.
##' @param ... additional arguments passed to \code{fn}.
##' @export
create_spline_matrix <- function(x, df=1, fn="ns",...){
    if (!is.numeric(df) || df<1) stop("'df' must be postive integer.")
    if (df==1){
        Hx <- matrix(x)
    } else {
        fn <- get(fn)
        Hx <- fn(x=x,
                 df=df,
                 ...)
    }
    Hx
}

#' @rdname create_spline_matrix
#' @export
create_spline <- function(x, df=1, fn="ns",...){
    spline <- create_spline_matrix(x, df=df, fn=fn, ...)
    out <- list(x=x,
                spline=spline,
                spline_mean=colMeans(spline))
    out
}


#' @rdname create_spline_matrix
#' @param xseq Sequence of values to use for creating time spline.
#' @param xrange Range (min and max) of values to use for creating spline. Only used if \code{xseq} not provided.
#' @param by Interval of generated sequence.
##' @export
##' @seealso \code{\link{sim_update_time_splines}}
create_spline_range <- function(df, xseq=NULL, xrange=NULL, by=1, fn="ns", ...){
    if(is.null(xseq)) {
        if (is.null(xrange)){
            stop("Either 'xseq' or 'xrange' required.")
        } else {
            xseq <- seq(xrange[1], to=xrange[2], by=by)
        }
    }
    out <- create_spline(xseq, df=df, fn=fn, ...)
    out
}

##' @rdname create_spline_matrix
##' @param spline Spline object from which to predict.
##' @param newx New values for which splines should be generated
##' @param center,scale Should spline values be centered and/or scaled? Default is \code{FALSE}. If a vector is provided, it is pased to \code{\link{scale}}.
##' @export
##' @importFrom utils getS3method
predict_spline <- function(spline, newx, center=FALSE, scale=FALSE){
    spline_class <- class(spline)
    if(spline_class[1]=="matrix" && length(spline_class)==1) {
        # NOT checking 'inherits' here, since most splines are matrices, but not vice versa
        # THen odf = 1
        if (ncol(spline)!=1) stop("Error: spline has >1 column but is class matrix only")
        H <- as.matrix(newx)
    } else {
        predfn <- utils::getS3method(f="predict", class(spline))
        H <- predfn(spline, newx=newx)
    }

    if(any(as.logical(center)) || any(as.logical(scale))){
        H <- scale(H, center=center,scale=scale)
    }
    H
}


##' @title Add time spline to a standata object
##' @description Updates a standata object with a matrix of time or exposure splines
##' @param standata the standata object to be updated
##' @param Ht,Hx the matrix of time or exposure splines to be added
##' @details \code{add_Ht_standata} updates the \code{Ht}, \code{timedf}, and \code{Ht_attr} elements of \code{standata}. \code{add_Hx_standata} updates the \code{Hx}, \code{xdf}, and \code{Hx_attr} elements of \code{standata}.
##' @seealso \code{\link{create_standata_exposure}}, \code{\link{create_standata_outcome}}, \code{\link{expsim_update_parameter}}, \code{\link{outsim_update_parameter}}, \code{\link{create_spline_range}}
##' @export
add_Ht_standata <- function(standata, Ht){
    if(!inherits(Ht, c("matrix","numeric", "integer"))) stop("'Ht' should be a matrix or vector.")
    standata$Ht <- as.matrix(Ht)
    standata$timedf <- ncol(standata$Ht)
    standata$Ht_attributes <- attributes(Ht)
    standata
}

#' @rdname add_Ht_standata
#' @export
add_Hx_standata <- function(standata, Hx){
    if(!inherits(Hx, c("matrix","numeric", "integer"))) stop("'Hx' should be a matrix or vector.")
    standata$Hx <- as.matrix(Hx)
    standata$xdf <- ncol(standata$Hx)
    standata$Hx_attributes <- attributes(Hx)
    standata
}

#' @title Add Duplicated Observations Flag
#' @description Adds a column indicating duplicated rows
#' @param data Data frame to check.
#' @param dupvars Character vector providing the names of variables for which duplicated observations should be checked.
#' @details This is a wrapper around \code{\link{duplicated}}.
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



