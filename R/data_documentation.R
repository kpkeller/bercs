#' Example data giving concentration, probabilities of respiratory infection, and study information
#'
#'
#' @format A data frame with 200 rows and 6 columns:
#' \describe{
#'   \item{conc}{concentration of exposure}
#'   \item{case}{whether or not an infection was observed}
#'   \item{unit_id}{The id corresponding to each individual in the study}
#'   \item{true_prob}{true probability of observing a case}
#'   \item{true_logit_prob}{the logit of true_prob}
#'   \item{study}{The study from which the data was taken}
#' }
#' @details detailssss here
'casedataA'

#' Example data giving concentration, probabilities of respiratory infection, and study information
#'
#' @rdname casedataA
#'
'casedataB'

#' Example data giving a linear relationship between conc and sbp
#'
#' @rdname casedataA
#' @format A data frame with 200 rows and 6 columns:
#' \describe{
#'   \item{sbp}{systolic blood pressure}
#'   ...
#' }
#' @examples
#' data(casedataC)
#' plot(x=casedataC$conc, y=casedataC$sbp)
'casedataC'

#' Example data giving a logarithmic relationship between conc and sbp
#'
#' @rdname casedataA
'casedataD'

#' Example data giving an exponential relationship between conc and sbp
#'
#' @rdname casedataA
'casedataE'

#' Example data giving a piecewise relationship between conc and sbp
#'
#' @rdname casedataA
'casedataF'





