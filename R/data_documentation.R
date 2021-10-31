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
'casedataA'

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
'casedataB'

#' Example data giving a linear relationship between conc and sbp
#'
#'
#' @format A data frame with 200 rows and 6 columns:
#' \describe{
#'   \item{conc}{concentration of exposure}
#'   \item{sbp}{systolic blood pressure}
#'   ...
#' }
#' @examples
#' data(casedataC)
#' plot(x=casedataC$conc, y=casedataC$sbp)
'casedataC'

#' Example data giving a logarithmic relationship between conc and sbp
#'
#'
#' @format A data frame with 200 rows and 6 columns:
#' \describe{
#'   \item{conc}{concentration of exposure}
#'   \item{sbp}{systolic blood pressure}
#'   ...
#' }
#' @examples
#' data(casedataD)
#' plot(x=casedataD$conc, y=casedataD$sbp)
'casedataD'

#' Example data giving an exponential relationship between conc and sbp
#'
#'
#' @format A data frame with 200 rows and 6 columns:
#' \describe{
#'   \item{conc}{concentration of exposure}
#'   \item{sbp}{systolic blood pressure}
#'   ...
#' }
#' @examples
#' data(casedataE)
#' plot(x=casedataE$conc, y=casedataE$sbp)
'casedataE'

#' Example data giving a piecewise relationship between conc and sbp
#'
#'
#' @format A data frame with 200 rows and 6 columns:
#' \describe{
#'   \item{conc}{concentration of exposure}
#'   \item{sbp}{systolic blood pressure}
#'   ...
#' }
#' @examples
#' data(casedataF)
#' plot(x=casedataF$conc, y=casedataF$sbp)
'casedataF'


