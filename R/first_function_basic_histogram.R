library(bercs)
library(ggplot2)
str(casedataA)
#casedataA



#' @title get_conc_hist
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
        geom_histogram(bins=20, aes(vector)) +
        ggtitle('Frequency of Concentrations') +
        labs(y='Frequency', x='Concentration')

}



