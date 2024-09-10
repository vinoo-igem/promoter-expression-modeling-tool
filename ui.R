library(shiny)
library(ggplot2)
library(bslib)

dataset <- diamonds

fluidPage(
  
  titlePanel("DBTL Bootcamp: Modeling Constiuitive Promoter Expression"),
  
  layout_columns(
    card( 
      card_header("Model"),
      withMathJax(),
      div(
        p("Transcription"),
        p("$$\\frac{d[mRNA]}{dt} = k_{1}{[Gene]} - d_{1}{[mRNA]}$$"),
        p("Translation"),
        p("$$\\frac{d[Protein]}{dt} = k_{2}{[mRNA]} - d_{2}{[Protein]}$$")
      )
    ),
    card( 
      card_header("Promoter Strength Data (k1)"),
      p("enter the strength data of your promoter"),
      fluidRow(
        column(8, textInput('pro1Name', 'Name/Label', value = 'J23100')),
        column(4, numericInput('pro1Str','strength',1.00,min = 0,max = NA,step = 0.01))
      ),
      fluidRow(
        column(8, textInput('pro2Name', 'Name/Label', value = 'J23101')),
        column(4, numericInput('pro2Str','strength',0.70,min = 0,max = NA,step = 0.01))
      ),
      fluidRow(
        column(8, textInput('pro3Name', 'Name/Label', value = 'J23106')),
        column(4, numericInput('pro3Str','strength',0.47,min = 0,max = NA,step = 0.01))
      ),
      fluidRow(
        column(8, textInput('pro4Name', 'Name/Label', value = 'J23116')),
        column(4, numericInput('pro4Str','strength',0.16,min = 0,max = NA,step = 0.01))
      ),
      fluidRow(
        column(8, textInput('pro5Name', 'Name/Label', value = 'J23117')),
        column(4, numericInput('pro5Str','strength',0.06,min = 0,max = NA,step = 0.01))
      )
    ),
    card( 
      card_header("Other Model Parameters & States"),
      p("these parameters and states will be used across all modeled constructs"),
      numericInput('mrnaDeg','mRNA degradation (d1)',0.23,min = 0,max = NA,step = 0.01,width = NULL),
      numericInput('rbsStr','RBS strength (k2)',8.23,min = 0,max = NA,step = 0.01,width = NULL),
      numericInput('protDeg','protein degradation (d2)',0.02,min = 0,max = NA,step = 0.01,width = NULL),
      numericInput('copyNumber','plasmid copy number (Gene)',17,min = 0,max = NA,step = 1,width = NULL)
    )
  ),
  mainPanel(
    plotOutput('plot')
  )
)