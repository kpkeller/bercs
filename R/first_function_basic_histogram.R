library(bercs)
library(ggplot2)
str(casedataA)



#' @title get_conc_hist
#'
#' @param vector a vector of concentrations from a study
#' @param ... multiple data frames
#'
#' @return gives a histogram of concentrations
#' @details This function gives you a historam for the counts of concentrations
#' @export
#'
#' @examples
get_conc_hist <- function(vector, ...){
    num_pars <- length(match.call())-1

    if(num_pars ==1){
    ggplot() +
        geom_histogram(bins=20, aes(vector)) +
        ggtitle('Frequency of Concentrations') +
        labs(y='Frequency', x='Concentration')
    }

    else{
        combined_data <- rbind(...)
        ggplot() +
            geom_histogram(bins=20, aes(combined_data$conc)) +
            ggtitle('Frequency of Concentrations') +
            labs(y='Frequency', x='Concentration')


    }

}





