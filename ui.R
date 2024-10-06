library(shiny)
library(ggplot2)
library(bslib)

dataset <- diamonds

fluidPage(
  
  titlePanel("DBTL Bootcamp: Modeling Constitutive Promoter Protein Expression"),
  p("The following are the ordinary differential equations (ODEs) we are 
        using in the iGEM DBTL Bootcamp to model protein expression under 
        constitutive promoters. For more information, please see the modeling section in the DBTL Bootcamp 
        documentation."),
  sidebarLayout(
    sidebarPanel(
      card( 
        card_header(tags$h3("Parameters & States")),
        p("these parameters and states are used across all modeled constructs"),
        numericInput('copyNumber','plasmid copy number (Gene)',17,min = 0,max = NA,step = 1,width = NULL),
        numericInput('mrnaDeg','mRNA degradation (d1)',0.23,min = 0,max = NA,step = 0.01,width = NULL),
        numericInput('rbsStr','RBS strength (k2)',8.23,min = 0,max = NA,step = 0.01,width = NULL),
        numericInput('protDeg','protein degradation (d2)',0.02,min = 0,max = NA,step = 0.01,width = NULL)
      ),
      card( 
        card_header(tags$h3("Promoter Strength Data (k1)")),
        p("enter the name and strength (str) data of your promoter"),
        fluidRow(
          column(8, textInput('pro1Name', 'Name/Label', value = 'J23100')),
          column(4, numericInput('pro1Str','str',1.00,min = 0,max = NA,step = 0.01))
        ),
        fluidRow(
          column(8, textInput('pro2Name', 'Name/Label', value = 'J23101')),
          column(4, numericInput('pro2Str','str',0.70,min = 0,max = NA,step = 0.01))
        ),
        fluidRow(
          column(8, textInput('pro3Name', 'Name/Label', value = 'J23106')),
          column(4, numericInput('pro3Str','str',0.47,min = 0,max = NA,step = 0.01))
        ),
        fluidRow(
          column(8, textInput('pro4Name', 'Name/Label', value = 'J23116')),
          column(4, numericInput('pro4Str','str',0.16,min = 0,max = NA,step = 0.01))
        ),
        fluidRow(
          column(8, textInput('pro5Name', 'Name/Label', value = 'J23117')),
          column(4, numericInput('pro5Str','str',0.06,min = 0,max = NA,step = 0.01))
        ),
        actionButton(
          inputId = "submitData",
          label = "Submit")
      )
    ),
    mainPanel(
      card( 
        card_header(tags$h3("Model")),
        fluidRow(
          column(4, img(src = "transcription-unit-k-d.png", style = "width: 100%; padding: 0;")),
          column(8,
                 withMathJax(),
                 div(
                   p("Transcription"),
                   p("$$\\frac{d[mRNA]}{dt} = k_{1}{[Gene]} - d_{1}{[mRNA]}$$"),
                   p("Translation"),
                   p("$$\\frac{d[Protein]}{dt} = k_{2}{[mRNA]} - d_{2}{[Protein]}$$")
                 )
          )
        )
      ),
      fluidRow(
        column(6, plotOutput('plot')),
        column(6, plotOutput('chart'))
      )


    )
  )
)