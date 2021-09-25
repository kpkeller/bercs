library(bercs)
library(ggplot2)
str(casedataA)
#casedataA



#' get_conc_hist
#'
#' @param vector
#'
#' @return gives a histogram of concentrations
#' @details This function gives you a historam for the counts of concentrations
#' @export
#'
#' @examples
get_conc_hist <- function(vector){
    ggplot() +
        geom_histogram(aes(vector)) +
        ggtitle('Frequency of Concentrations') +
        labs(y='Frequency', x='Concentration')

}

get_conc_hist(casedataA$conc)
?get_conc_hist
