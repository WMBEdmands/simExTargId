#' xcms results peak monitor and shiny application
#' @param xcmsOutput character or data.frame. Full path to xcms results file (or an XCMS output file peak table in the form of a data.frame) if argument is not supplied a tcltk window will open.
#' @param obsNames character vector of sample/QC/Blank names to identify these columns in the xcms output table. The obsNames argument should be pre-sorted in to the acquisition run order.
#' @param qcNames character vector of analytical replicate local QC names to identify these columns in the xcms output table. The qcNames argument should be pre-sorted in to the acquisition run order. These identificatory sample names are critical for CV\% calculation for example.
#' @param metab character or data.frame. Full path to a csv file of a list of metabolites you wish to monitor (or a data.frame) if argument is not supplied a tcltk window will open.
#' This list will be used to identify peak groups in the xcms output table within
#' user defined ppm mass accuracy and retention time deviation parameters.
#' The table must consist of the following 3 columns minimum
#'
#' \enumerate{
#' \item name -- unique name of compound
#' \item mzmed -- expected mass-to-charge ratio
#' \item rtmed -- expected retention time (in seconds)
#' }
#'
#' @param ppm numeric ppm mass accuracy for matching monitored metabolites
#' to the xcms output table (default = 10 ppm).
#' @param rtdev numeric maximum deviation in seconds for matching monitored
#' metabolite list and the xcms output table peak groups (default = 30 seconds).
#' @param maxSignalAtt numeric percentage (i.e. values between 0-100\%, default=20\%)
#' maximum percentage last quality control signal attenuation from the quality controls median value.
#' If a monitored metabolite drops below this value then (i.e. -20\%) then it will be
#' flagged in the returned table.
#' @param maxBelow numeric (values between 0-100\%, default = 20 \%) maximum percentage of monitored
#' metabolites below the maximum signal attenutation argument (see argument maxSignalAtt). E.g. if 20\% of monitored metabolites contained in expected metabolite table have dropped
#' below 20\% median QC signal signal attenuation.
#' @param obsTypeStr character vector equal in length to obsNames to identify different run types
#' in the obsNames argument.
#' @param obsTypeCol character vector colours equal in length to obsNames to display for run types in application.
#' @param launchApp logical default is true. If true peakMonitor shiny app will open in your default web browser.
#' @param analysisDir the full path of the peakMonitor directory (output/peakMonitor)
#' containing the output files (.csv) from \code{\link{peakMonitor}}. This
#' argument is only required if the argument launchApp is TRUE. If the argument
#' is not supplied but required then a tcltk window will open and the directory
#' can be selected.
#' @return a list containing two elements named: \enumerate{
#' \item "monitMetab" a data.frame containing the metabolite peak monitoring results.
#' \item "percBelow" a logical whether the user parameterized maximum percentage of
#' monitored metabolites were below the maximum signal attenuation value. This value
#' is used by simExTargId to notify the user of any significant signal attenuation drop during
#' a run for example.
#' }
#' @export
peakMonitor <- function(xcmsOutput=NULL, obsNames=NULL, metab=NULL, qcNames=NULL,
                        ppm=10, rtdev=30, maxSignalAtt=20, percBelow=20,
                        obsTypeStr=rep('unknown_obsType', length(obsNames)),
                        obsTypeCol=rep('darkgrey', length(obsNames)),
                        launchApp=TRUE, analysisDir=NULL){

  if(launchApp == FALSE){
  # error handling
  if(is.null(obsNames)){
    stop('The names of the observation (i.e. samples, QCs, blanks) columns must be supplied')
  }
  if(is.null(qcNames)){
    stop('The names of the pooled quality control columns must be supplied')
  }

  if(is.null(xcmsOutput)){
    xcmsOutput <- tcltk::tclvalue(tcltk::tkgetOpenFile(title="select your xcms results file"))
  }

  if(is.character(xcmsOutput)){
  xcmsOutFileName <- basename(xcmsOutput)
  if(!grepl('\\.csv$', xcmsOutput)){
    stop('xcmsOutput file must be a .csv')
  }
    xcmsOutput <- as.data.frame(data.table::fread(xcmsOutput, header=TRUE, stringsAsFactors = FALSE))
  } else {
    xcmsOutFileName <- 'unknown'
  }

  if(!is.data.frame(xcmsOutput)){
    stop("xcmsOutput must be supplied as either a full path to an xcms output file or as a data.frame")
  }

  if(is.null(xcmsOutput$mzmed)){
    stop('xcmsOutput does not contain a column named "mzmed"')
  }
  if(is.null(xcmsOutput$rtmed)){
    stop('xcmsOutput does not contain a column named "rtmed"')
  }
  if(is.null(xcmsOutput$name)){
    gnames <- paste0('M', round(xcmsOutput$mzmed, 0), 'T',
                              round(xcmsOutput$rtmed, 0))
   if(any(dup <- duplicated(gnames))){
      for(dupname in unique(gnames[dup])) {
        dupidx <- which(gnames == dupname)
        gnames[dupidx] <- paste(gnames[dupidx], seq(along = dupidx), sep = "_")
     }
   }
   xcmsOutput$name <- gnames
  }

  if(is.null(metab)){
    metab <- tcltk::tclvalue(tcltk::tkgetOpenFile(title="select your metabolites to monitor file"))
  }
  if(is.character(metab)){
    if(!grepl('\\.csv$', metab)){
      stop('metabolite table file must be a .csv')
    }
    metab <- as.data.frame(data.table::fread(metab, header=TRUE, stringsAsFactors = FALSE))
  }

  if(!is.data.frame(metab)){
    stop("metabolite table must be supplied as either a full path to a csv file or as a data.frame")
  }
  reqCols <- c('name', 'mzmed', 'rtmed') %in% names(metab)
  if(any(reqCols == FALSE)){
    stop('The metabolite table must consist of a minimum of 3 column names:\n1. "name"\n2. "mzmed"\n3. "rtmed"\n')
  }
  # all unique compound names
  if(any(duplicated(metab$name))){
    stop('metabolite table names must be unique. If you have the same compound name twice (i.e. with different retention times place an additional unique integer at the end of the name. e.g. lysoPC (14:0) 1, lysoPC (14:0) 1')
  }
  if(length(obsTypeStr) != length(obsNames)){
    stop('argument obsTypeStr is not equal to obsName argument')
  }
  if(length(obsTypeCol) != length(obsNames)){
    stop('argument obsTypeCol is not equal to obsName argument')
  }
  # match metabs to xcmsOutput
  matchIdx <- as.numeric()
  metIdx <- as.numeric()

  for(i in 1:nrow(metab)){
  ppmDiff <- {{xcmsOutput$mzmed - metab$mzmed[i]}/xcmsOutput$mzmed} * 1E06
  mIdx <- abs(ppmDiff) <= ppm
  rIdx <- abs(xcmsOutput$rtmed - metab$rtmed[i]) <= rtdev
  mrIdx <- mIdx & rIdx
   if(any(mrIdx)){
   metIdx  <- c(metIdx, rep(i, sum(mrIdx)))
   matchIdx <- c(matchIdx, which(mrIdx))
   } else {
     warning('metabolite ', i, ' not found in xcmsOutput table (name = ', metab$name[i], ').\n', immediate. = TRUE)
   }
  }
  # expected metabolites
  exMetab <- cbind(metNo=as.integer(metIdx), metab[metIdx, , drop=FALSE])
  # remove : and /
  exMetab$name <- gsub('\\:|/', '.', exMetab$name)

  names(exMetab)[names(exMetab) %in% 'mzmed'] <- 'exp_mzmed'
  names(exMetab)[names(exMetab) %in% 'rtmed'] <- 'exp_rtmed'
  # observed metabolites
  obsMetab <- xcmsOutput[matchIdx, , drop=FALSE]
  exMetab$obs_mzmed <- obsMetab$mzmed
  exMetab$obs_rtmed <- obsMetab$rtmed

  exMetab$rtdev <- exMetab$exp_rtmed - obsMetab$rtmed
  exMetab$ppmDiff <- {{exMetab$exp_mzmed - obsMetab$mzmed}/exMetab$exp_mzmed} * 1E06

  # summary stats
  # cv%
  exMetab$qcCv <- apply(obsMetab[, qcNames, drop=FALSE], 1, function(x) {sd(x)/mean(x)} * 100)
  # mean/median intensity
  exMetab$qcMeanInt <- rowMeans(obsMetab[, qcNames, drop=FALSE])
  exMetab$qcMedianInt <- apply(obsMetab[, qcNames, drop=FALSE], 1, median)
  # curr mean/deviation
  exMetab$lastQcMeanDev <- round({{obsMetab[, qcNames[length(qcNames)]]/exMetab$qcMeanInt} - 1} * 100, 2)
  exMetab$lastQcMedianDev <- round({{obsMetab[, qcNames[length(qcNames)]]/exMetab$qcMedianInt} - 1} * 100, 2)
  exMetab$belowMeanSigAtt <- exMetab$lastQcMeanDev <= -maxSignalAtt
  exMetab$belowMedianSigAtt <- exMetab$lastQcMedianDev <= -maxSignalAtt
  # add sample data to exMetab table
  exMetab <- cbind(exMetab, obsMetab[, obsNames, drop=FALSE])
  # summary stats
  exMetab <- exMetab[order(exMetab$qcCv), , drop=FALSE]
  dupIdx <- duplicated(exMetab$metNo) == FALSE
  meanCv <- mean(exMetab$qcCv[dupIdx])
  sumMeanPeakArea <- sum(exMetab$qcMeanInt[dupIdx])
  sumMedianPeakArea <- sum(exMetab$qcMedianInt[dupIdx])
  lastQCMeanDev <- mean(exMetab$lastQcMeanDev[dupIdx])
  lastQCMedianDev <- mean(exMetab$lastQcMedianDev[dupIdx])

  # reorder exMetab
  exMetab <- exMetab[order(exMetab$metNo), , drop=FALSE]
  percBelowSigAtt <- {sum(exMetab$belowMedianSigAtt[dupIdx])/sum(dupIdx)} * 100
  percBelowLog <- percBelowSigAtt >= percBelow
  warnMessage <- ''
  if(percBelowLog){
    warnMessage <- paste0('Greater than ', percBelow,
                          '% of the monitored metabolite peak ',
                          'areas are currently attenuated by more than ',
                          maxSignalAtt,
                          '%.\n\nBased on the deviation of the last QC sample (',
                          qcNames[length(qcNames)],
                          ') from the median value of all previously collected QC samples (n=',
                          length(qcNames) - 1,').\n\nPlease check your experiment.\n')
    warning(warnMessage, immediate. = TRUE)
  }

  # make table summary
  summaryTable <-  c('Current mean CV(%): ', round(meanCv, 1), 'sum mean peak area: ',
                     prettyNum(round(sumMeanPeakArea, 1), big.mark = ','), 'sum median peak area: ',
                     prettyNum(round(sumMedianPeakArea, 1), big.mark = ','),
                     'last QC mean deviation (%): ', round(lastQCMeanDev, 1),
                     'last QC median deviation (%): ', round(lastQCMedianDev, 1),
                     paste0('percentage metabolite(s) below ', round(maxSignalAtt, 2),
                            '% peak attenuation'),
                     round(percBelowSigAtt, 2),
                     'number metabolites(s) detected:', length(unique(exMetab$metNo)),
                     'number metabolites in database (.csv):', nrow(metab))

  cat(paste0(summaryTable, collapse = '\n'))
  summaryTable <- matrix(summaryTable, byrow = TRUE, ncol=2)
  colnames(summaryTable) <- c('metric', 'result')
  # round necessary columns
  fourDigits <- c('exp_mzmed', 'obs_mzmed')
  twoDigits <- c('exp_rtmed', 'obs_rtmed', 'ppmDiff', 'rtdev', 'qcCv', 'qcMeanInt', 'qcMedianInt')
  exMetab[, fourDigits] <- apply(exMetab[, fourDigits], 2, function(x) round(x, 4))
  exMetab[, twoDigits] <- apply(exMetab[, twoDigits], 2, function(x) round(x, 2))
  # add colour and label info to table
  colLabs <- data.frame(matrix('', nrow=2, ncol=ncol(exMetab)), stringsAsFactors = FALSE)
  colnames(colLabs) <- colnames(exMetab)
  colLabs[1, obsNames] <- obsTypeStr
  colLabs[2, obsNames] <- obsTypeCol
  exMetab <- rbind(colLabs, exMetab)

  return(list(monitMetab=exMetab, summaryTable=summaryTable,
              warnMessage=warnMessage))

  # if req launch shiny app
  } else if(launchApp){

    if(is.null(analysisDir)){
      message("window opened to select analysis (/peakMonitor) results directory containing the results (.csv) files...")
      flush.console()

      analysisDir <- tcltk::tk_choose.dir(default = "",
                                          caption = "Select the analysis (/peakMonitor) results directory containing the results (.csv) files...")
    }
    analysisDir <<- analysisDir
    if(!require(shiny)){
      stop('The package shiny must be installed to use the peakMonitor application...')
    }

    appDir <<- system.file("shiny-apps", "peakMonitor", package = "simExTargId")
    if (appDir == "") {
      stop("Could not find example directory. Try re-installing `simExTargId`.", call. = FALSE)
    }

   # calculate mean all compounds
   source(paste0(appDir, '/global.R'))
   source(paste0(appDir, '/ui.R'))
   source(paste0(appDir, '/server.R'))
   shiny::runApp(list(ui = peakMonitUi, server = peakMonitServer),
                 launch.browser = TRUE)
  }
} # end function
