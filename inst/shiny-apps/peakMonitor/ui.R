# UI
peakMonitUi <- shinyUI(fluidPage(
  titlePanel(paste0("Peak Monitor app (", basename(xcmsOutFileName), ")")),
  sidebarLayout(
    sidebarPanel(h4("Peak Monitor options:"),
                 
                 selectInput("Compound",label=tags$b("Select the name of the compound to display :"),choices=c("MeanAllCompounds", exMetab$name)),
                 # selectInput("Result_type",label=tags$b("Select the xcms result type to display :"),choices=c(xcms_results_colnames)),
                 selectInput("Exclude",label=tags$b("Select any samples you wish to exclude:"), choices=c(obsNames), multiple=T)),
    
    mainPanel(h4("Peak monitoring table and plot"), width=8,
              # ggvis::ggvisOutput("exMetabPlot"),
              DT::dataTableOutput('exMetabTable'),
              shiny::plotOutput('exMetabPlot', width = "100%",  height = "600px",
                                click = "exMetabPlot_click",
                                brush = brushOpts(id = "exMetabPlot_brush")),
              DT::dataTableOutput('exMetabDataPoints'),
              # h5(strong("CV%")),
              # (textOutput("QC_CV_text")),
              # (textOutput("Sample_CV_text")),
              # h5(strong("Average")),
              # (textOutput("Av_text")),
              # (textOutput("Sample_Av_text")),
              # h5(strong("Mass accuracy")),
              # (textOutput("Mass_accuracy")),
              # h5(strong("Retention time")),
              # (textOutput("Retention_time")),
              # tags$br(),
              shiny::h4("PCA monitored metabolites"),
              shiny::radioButtons("pcaPlotType", label="PCA plot type:", 
                                  choices=c("scores", "loadings")),
              shiny::plotOutput('pcaPlot1', brush='pcaPlot1_brush',
                                width = "100%",  height = "600px")
              # showOutput("pcaPlot1","polycharts"),
              # showOutput("pcaPlot2","polycharts"),
              # h5(strong("PCA all features (loadings)")),
              # showOutput("PCA_plot_loadings","polycharts")
    )
    
  )
  
)) # end ui