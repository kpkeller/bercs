library(bercs)
library(ggplot2)
library(dplyr)

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
get_conc_hist <- function(datalist){


    for(i in 1:length(datalist)){
        if(i==1){
            datalist[[i]] %>%
                mutate(number = i) -> total_data
        }
        else{
            datalist[[i]] %>%
                mutate(number = i) -> new_data
            total_data <- rbind(total_data, new_data)
        }

    }

    total_data$number <- as.factor(total_data$number)
    ggplot(total_data) +
        geom_histogram(aes(conc, fill=study), alpha=0.5) +
        ggtitle('Frequency of Concentrations') +
        labs(y='Frequency', x='Concentration')



    # if(num_pars ==1){
    # ggplot() +
    #     geom_histogram(bins=20, aes(vector)) +
    #     ggtitle('Frequency of Concentrations') +
    #     labs(y='Frequency', x='Concentration')
    # }
    #
    # else{
    #     combined_data <- rbind(...)
    #     ggplot() +
    #         geom_histogram(bins=20, aes(combined_data$conc)) +
    #         ggtitle('Frequency of Concentrations') +
    #         labs(y='Frequency', x='Concentration')
    #
    #
    # }

}




