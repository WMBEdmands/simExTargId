#' Simultaneous MS1 profiling and target identification visualization
#'
#' @description Visualization of deconvoluted data obtained during an MS1-level
#' profiling experiment. All statResults (.csv) files within an analysis directory
#' will be read and can be visualized as a volcano plot and as reactive
#' tables indicating the direction of the median fold change. This allows the
#' experimenter to visualize the emergence of statistically significant
#' pseudospectra features before the completion of a MS1-level profiling experiment.
#' minimum fold change and maximum multiple testing adjusted p-/q-values (\code{\link{coVarStatType}}) can be
#' set by the experimenter. When a list of possible targets has been determined
#' the table can be downloaded/ saved as a Csv file. This list of statistically
#' significant targets can therefore be readily included in a MS/MS fragmentation
#' experiment before the end of a profiling experiments completion.
#'
#' @param analysisDir (i.e. 04.stats)the full path of the analysis directory containing the
#'  (.csv) files output from \code{\link{simExTargId}}.
#'
#' @seealso \code{\link{simExTargId}}, \code{\link{rtCorrClust}}.
#'
#' @export
targIdShiny <- function(analysisDir=NULL){

  if(is.null(analysisDir)){
    message("window opened to select analysis (/04.stats) results directory containing the results (.csv) files...")
    flush.console()

    analysisDir <- tcltk::tk_choose.dir(default = "",
                                        caption = "Select the analysis (/04.stats) results directory containing the results (.csv) files...")
  }

  # identify results files
  resFileNames <- list.files(analysisDir, full.names=T, pattern="\\.csv",
                              recursive=T)

  # add all results files into list
  message("Reading ", length(resFileNames), " (.csv) results files please wait..")
  flush.console()
  resFiles <- lapply(resFileNames, read.csv, header=T, stringsAsFactors=F)
  resFileNames <- gsub("statResults_|\\.csv$", "", basename(resFileNames))
  names(resFiles) <- resFileNames
  # id p value columns
  FCPnames <- lapply(resFiles, function(x){
    pValNamesIndx <- grep("^p\\.value", colnames(x))
    pValNames <- colnames(x)[pValNamesIndx]
    pAdjNamesIndx <- grep("^Adj\\.p\\.value", colnames(x))
    pAdjNames <- colnames(x)[pAdjNamesIndx]
    FCnamesIndx <- grep("FoldChange", colnames(x))
    FCnames <- colnames(x)[FCnamesIndx]
    list(pValNamesIndx=pValNamesIndx, pValNames=pValNames,
         pAdjNamesIndx=pAdjNamesIndx, pAdjNames=pAdjNames,
         FCnamesIndx=FCnamesIndx, FCnames=FCnames)})

  # covariate names
  coVarNames <- unique(gsub("^p\\.value_", "", FCPnames[[1]]$pValNames))
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
                                                      choices=c("marginal p-value", "multiple testing adjusted")),
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
  targIdserver <- shiny::shinyServer(function(input, output, session){

    observe({ if(!is.null(input$Results_file)){
    resFileIndx <- which(names(resFiles) == input$Results_file)
    resFile <- resFiles[[resFileIndx]]
    FCPnames.tmp <- FCPnames[[resFileIndx]]
      if(input$multTest == "marginal p-value"){
        pvalIndx <- FCPnames.tmp$pValNames[grep(input$co_variate, FCPnames.tmp$pValNames)]
        log2FC <- log2(resFile[, FCPnames.tmp$FCnames[grep(input$co_variate, FCPnames.tmp$FCnames)]])
      } else {
        pvalIndx <- FCPnames.tmp$pAdjNames[grep(input$co_variate, FCPnames.tmp$pAdjNames)]
        log2FC <- log2(resFile[, FCPnames.tmp$FCnames[grep(input$co_variate, FCPnames.tmp$FCnames)]])
      }

    pvals_tmp <- resFile[, pvalIndx]
   # names(log2FC) <- resFile$RtCorrClust
    names(pvals_tmp) <- resFile$RtCorrClust

    # identify negative Fc and p-value
   if(!is.null(dim(log2FC))){
     # if anova results then if any above min Foldchange and take
     negFClogi <- log2FC < (-log2(input$minFC))
     posFClogi <- log2FC > log2(input$minFC)
     negFClogi <- rowSums(negFClogi) >= 1
     posFClogi <- rowSums(posFClogi) >= 1
     log2FC <- apply(log2FC, 1, max)
   } else {
    negFClogi <- log2FC < (-log2(input$minFC))
    posFClogi <- log2FC > log2(input$minFC)
   }
    pValLogi <- pvals_tmp < input$minPval
    # column display index
    colDispIndx <- grep("Feat|RtCorrClust|meanInt|mostIntSamp_|bestSampMSMS",
                        colnames(resFile))
    FCpIndx.tmp <- c(FCPnames.tmp$pValNamesIndx, FCPnames.tmp$pAdjNamesIndx,
                     FCPnames.tmp$FCnamesIndx)
    coVarIndx.tmp <- FCpIndx.tmp[grep(input$co_variate, colnames(resFile)[FCpIndx.tmp])]

    colDispIndx <- c(1:(min(FCpIndx.tmp) - 1), coVarIndx.tmp, colDispIndx)
    # volcano plot
    output$volcanoPlot <- shiny::renderPlot({
            MetMSLine::volcanoPlot(folds=log2FC, pvals=pvals_tmp,
                                   plimit=input$minPval, fclimit=input$minFC,
                                   cexcutoff = 1.5, cexlab = 1.5,
                                   ylab="-log10 p-value",
                                   xlab="log2 fold change")
    })
    # negFC table
    output$negFCtable <- DT::renderDataTable({
                       negFCtable.tmp <- resFile[negFClogi & pValLogi, colDispIndx]
                       return(negFCtable.tmp)
                       }, options = list(pageLength = 10))
    # posFC table
    output$posFCtable <- DT::renderDataTable({
                       posFCtable.tmp <- resFile[posFClogi & pValLogi, colDispIndx]
                       return(posFCtable.tmp)
                       }, options = list(pageLength = 10))

    # nSignifFeatures table
    output$nSignifFeat <- shiny::renderTable({
      sumNegFC <- sum(negFClogi & pValLogi)
      sumPosFC <- sum(posFClogi & pValLogi)
      nSignifFeat.df <- data.frame(sumNegFC, sumPosFC, stringsAsFactors=F)
      colnames(nSignifFeat.df) <- paste0("n significant features ", c("NEG FC","POS FC"))
      return(nSignifFeat.df)
    })
    # download data
    output$downloadData <- shiny::downloadHandler(
      filename = function() { paste0("TargId_",input$Results_file, "_",
                                     input$co_variate,
                                     ifelse(input$multTest == "marginal p-value",
                                            "", "_multTestAdjusted"), '.zip') },
      content = function(file) {
        Targets <- rbind(resFile[negFClogi & pValLogi, colDispIndx],
                         resFile[posFClogi & pValLogi, colDispIndx])
        tmpdir <- tempdir()
        setwd(tempdir())
        print(tempdir())

        csvIndiv <- by(Targets, Targets$bestSampMSMS, function(x){
          write.csv(x, file=paste0("TargId_",input$Results_file, "_",
                                   input$co_variate,
                                   ifelse(input$multTest == "marginal p-value",
                                          "", "_multTestAdjusted"), "_",
                                   unique(x[, "bestSampMSMS"]), ".csv"),
                    row.names=F)})

        fs <- paste0("TargId_",input$Results_file, "_",
                     input$co_variate,
                     ifelse(input$multTest == "marginal p-value",
                            "", "_multTestAdjusted"), "_",
                     unique(Targets[, "bestSampMSMS"]), ".csv")
        print(fs)

        zip(zipfile=file, files=fs)
      },
      contentType = "application/zip")
    }}) # end observer 1

  }) # END compMS2server

  shiny::runApp(list(ui=targIdui, server=targIdserver), launch.browser=T)

} # end function
