library(shiny)
library(ggplot2)
library(reshape2)
library(deSolve)

function(input, output) {

  observeEvent( input$submitData, {
    
    # NOTE: need to become dynamic, with user selectable # of inputs
    # collect the promoter strengths
    promoter_strengths <- c(
      input$pro1Str,
      input$pro2Str,
      input$pro3Str,
      input$pro4Str,
      input$pro5Str
    )
    
    # collect the names of the promoters
    promoter_names <- c(
      input$pro1Name,
      input$pro2Name,
      input$pro3Name,
      input$pro4Name,
      input$pro5Name
    )
    
    # NOTE: needs better error handling, simple for now
    # grab indices where promoter name is not empty, and strength exists and is greater than 0
    valid_indices <- which(promoter_names != "" & !is.na(promoter_strengths) & promoter_strengths > 0)
    
    # only keep valid promoter names and strengths
    valid_promoter_strengths <- promoter_strengths[valid_indices]
    valid_promoter_names <- promoter_names[valid_indices]
    
    # create promoter name and strength as key-value pairs
    # NOTE: probably an easier way to do this, coupled with error handling and above code
    k1_values <- setNames(valid_promoter_strengths, valid_promoter_names)
    
    # all default parameters and states taken from https://youtu.be/K0P1KVk_hDo?si=5K72owA1k8SNi9u7&t=1046
    # NOTE: these will need to be checked
    parameters <- c(
      k1 = NA, # initialize, but values will be added from k1_values in the loop
      d1 = input$mrnaDeg,
      k2 = input$rbsStr,
      d2 = input$protDeg
    )
    
    state <- c(gene = input$copyNumber, mRNA = 0, protein = 0)
    
    all_results_df <- data.frame()
    
    # loop through each promoter strength value
    for (name in names(k1_values)) {
      
      k1_value <- k1_values[[name]]  # access reactive k1 values
      
      # and assign promoter strength value to parameters vector
      parameters["k1"] <- k1_value
      
      derivative <- function(times, state, parameters) {
        with(as.list(c(state, parameters)), {
          # Rate of change
          dGene <- 0
          dmRNA <- k1 * gene - d1 * mRNA
          dProtein <- k2 * mRNA - d2 * protein
          return(list(c(dGene, dmRNA, dProtein)))
        })
      }
      
      times <- seq(0, 180, by = 0.1) # time points in minutes
      
      # solve the ODE
      out <- ode(y = state, times = times, func = derivative, parms = parameters)
      
      # convert ODE output to a data frame
      out_df <- as.data.frame(out)
      
      # add column and name to identify the k1 (promoter)
      out_df$k1 <- name
      
      # add the results to the data frame
      all_results_df <- rbind(all_results_df, out_df)
    }
    
    # melt the data frame for ggplot2
    results_melt <- melt(all_results_df, id.vars = c("time", "k1"))
    
    # filter out non-positive values, logarithms of zero or negative numbers are undefined
    results_melt <- results_melt[results_melt$value > 0, ]
    
    output$plot <- renderPlot({
      ggplot(results_melt, aes(x = time, y = value, color = variable, linetype = k1)) +
        geom_line() +
        scale_y_log10() +
        labs(title = "Gene Expression by Promoter Strength",
             x = "Time (min)",
             y = "Concentration (log scale)",
             color = "Component",
             linetype = "Promoter") +
        theme_minimal()
      
    }, height = 500)
    
    # pull only protein and max time-point data for chart
    bar_chart_data <- results_melt[results_melt$variable == "protein" & results_melt$time == max(results_melt$time), ]
    
    output$chart <- renderPlot({
      ggplot(bar_chart_data, aes(x = k1, y = value, fill = k1)) +
        geom_bar(stat = "identity") +
        labs(title = "Total Protein Concentration by Promoter",
             x = "Promoter",
             y = "Total Protein Concentration") +
        theme_minimal() +
        theme(legend.position = "none")
      
    }, height = 500)
    
    }, ignoreNULL = FALSE)
  
}
