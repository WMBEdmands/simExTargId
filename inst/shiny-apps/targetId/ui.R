# Shiny UI    d
targIdui <- shiny::shinyUI(shiny::fluidPage(
  shiny::titlePanel("Simultaneous Experiment - MS/MS Target Identification "),
  shiny::fluidRow(shiny::column(shiny::h4("Options"),width=3,
                                shiny::selectInput("Results_file", label=shiny::tags$b("select the results output file:"),
                                                   names(resFiles)),
                                shiny::selectInput("co_variate", label=shiny::tags$b("select co-variate:"),
                                                   coVarNames),
                                shiny::numericInput("minFC",label=shiny::tags$b("minimum fold change : "),
                                                    value=2, min=1, max=20, step=0.1),
                                shiny::numericInput("minPval",label=shiny::tags$b("maximum p-value : "),
                                                    value=0.05, min=0.001, max=1, step=0.001),
                                shiny::radioButtons("multTest", label=shiny::tags$b("p value multiple testing adjustment : "),
                                                    choices=c("raw p-value", "multiple testing adjusted")),
                                shiny::tableOutput("nSignifFeat"),
                                shiny::downloadButton('downloadData', 'Download zip file of filtered target tables (.csv)')
  ), # end column 1
  shiny::column(width=8,
                shiny::tabsetPanel(
                  shiny::tabPanel("volcano plot",
                                  shiny::plotOutput("volcanoPlot", height="600px")),
                  shiny::tabPanel("Targets negative fold change", DT::dataTableOutput(outputId="negFCtable")),
                  shiny::tabPanel("Targets positive fold change", DT::dataTableOutput(outputId="posFCtable"))
                ) # end tabset panel
  ) # end column 2
  ) # end fluid Row 1
) # end fluidPage
) # end CompMS2 shiny UI
