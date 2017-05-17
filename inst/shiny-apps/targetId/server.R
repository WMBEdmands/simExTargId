
targIdserver <- shiny::shinyServer(function(input, output, session){

  session$onSessionEnded(function(){
    shiny::stopApp()})

  observe({ if(!is.null(input$Results_file)){
    resFileIndx <- which(names(resFiles) == input$Results_file)
    resFile <- resFiles[[resFileIndx]]
    FCPnames.tmp <- FCPnames[[resFileIndx]]
    if(input$multTest == "raw p-value"){
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
                                     ifelse(input$multTest == "raw p-value",
                                            "", "_multTestAdjusted"), '.zip') },
      content = function(file) {
        Targets <- rbind(resFile[negFClogi & pValLogi, colDispIndx],
                         resFile[posFClogi & pValLogi, colDispIndx])
        tmpdir <- tempdir()
        setwd(tempdir())
        print(tempdir())

        # csvIndiv <- by(Targets, Targets$bestSampMSMS, function(x){
        #   write.csv(x, file=paste0("TargId_",input$Results_file, "_",
        #                            input$co_variate,
        #                            ifelse(input$multTest == "raw p-value",
        #                                   "", "_multTestAdjusted"), "_",
        #                            unique(x[, "bestSampMSMS"]), ".csv"),
        #             row.names=FALSE)})
        write.csv(Targets, paste0("TargId_",input$Results_file, "_",
                                  input$co_variate,
                                  ifelse(input$multTest == "raw p-value",
                                         "", "_multTestAdjusted"), '.csv'),
                  row.names = FALSE)
        fs <- paste0("TargId_",input$Results_file, "_",
                     input$co_variate,
                     ifelse(input$multTest == "raw p-value",
                            "", "_multTestAdjusted"), '.csv') #"_",
                     #unique(Targets[, "bestSampMSMS"]), ".csv")
        print(fs)

        zip(zipfile=file, files=fs)
      },
      contentType = "application/zip")
  }}) # end observer 1

}) # END compMS2server
