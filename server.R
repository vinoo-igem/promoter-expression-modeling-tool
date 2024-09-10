library(shiny)
library(ggplot2)
library(reshape2)
library(deSolve)

function(input, output) {

  k1_values1 <- reactive({
    c(
      pro1Name = input$pro1Str,
      pro2Name = input$pro2Str,
      pro3Name = input$pro3Str,
      pro4Name = input$pro4Str,
      pro5Name = input$pro5Str
    )
  })
  
  
  # Collect the k1 values into a list
  promoter_strengths <- reactive({list(
    input$pro1Str,
    input$pro2Str,
    input$pro3Str,
    input$pro4Str,
    input$pro5Str
  )})
  
  # Collect the names of the promoters
  promoter_names <- reactive({c(
    input$pro1Name,
    input$pro2Name,
    input$pro3Name,
    input$pro4Name,
    input$pro5Name
  )})
  
  k1_values <- reactive(setNames(promoter_strengths(), promoter_names()))
  
  
  # reactive expression for solving the ODE and creating the data frame
  all_results <- reactive({
    all_results_df <- data.frame()
    
    # parameters and states taken from https://youtu.be/K0P1KVk_hDo?si=5K72owA1k8SNi9u7&t=1046
    # these will need to be checked
    parameters <- c(
      k1 = NA, # initialize, but values will be added from k1_values in the loop
      d1 = input$mrnaDeg,
      k2 = input$rbsStr,
      d2 = input$protDeg
    )
    
    state <- c(gene = input$copyNumber, mRNA = 0, protein = 0)
    
    # loop through each promoter strength value
    for (name in names(k1_values())) {
      
      k1_value <- k1_values()[[name]]  # access reactive k1 values
      
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
      
      # Ensure time column exists
      if (!"time" %in% colnames(out_df)) {
        out_df$time <- times  # Add time column if missing
      }
      
      # add column and name to identify the k1 (promoter)
      out_df$k1 <- name
      
      # add the results to our data frame
      all_results_df <- rbind(all_results_df, out_df)
    }
    
    return(all_results_df)
  })
  
  output$plot <- renderPlot({
    # Get all results
    results <- all_results()
    
    # melt the data frame for ggplot2
    results_melt <- melt(results, id.vars = c("time", "k1"))
    
    # filter out non-positive values, logarithms of zero or negative numbers are undefined
    results_melt <- results_melt[results_melt$value > 0, ]
    
    ggplot(results_melt, aes(x = time, y = value, color = variable, linetype = k1)) +
      geom_line() +
      scale_y_log10() +
      labs(title = "Gene Expression by Promoter Strength",
           x = "Time (min)",
           y = "Concentration (log scale)",
           color = "Component",
           linetype = "Promoter") +
      theme_minimal()
    
  }, height = 700)
  
}
