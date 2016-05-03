#' simExTargId
#' 
#' @description This function seeks to remove the gap between metabolomic 
#' MS1 profiling experiments and discovery of statistically relevant targets for 
#' MSn fragmentation based identification.  
#' 
#' @param rawDir character full path name of raw data directory into which raw data will 
#' eventually/ has already been started to be written. An "_analysis" directory
#' will be automatically created into which all converted files (.mzXML) and
#' results output will be saved. Fifteen minutes after the last modification of a new
#' raw data file a conversion to the mzXML format and initial peak-picking will
#' occur (\code{\link{xcmsSet}}). 
#' @param analysisDir full path name of analysis directory, all mzXML, peak-picking
#' and results files will be saved here.
#' @param coVar covariate table can be supplied as either a full path name to a 
#' .csv file (see details) or as a \code{\link{data.frame}}. 
#' @param nSlaves numeric number of computer cores for parallel computation.
#' @param minFiles the minimum number of raw data file converted to mzXML files
#' before commencing subsequent steps of XCMS processing, pre-processing, PCA and
#' statistical analysis. Default = 10. 
#' @param centroid do the raw data files need to be centroided during conversion by MSConvert
#' \url{http://proteowizard.sourceforge.net/downloads.shtml}. NB. centroiding of data
#' is necessary for utilization of the "centWave" algorithm of xcms (\code{\link{findPeaks.centWave-methods}}).
#' @param MS2 logical should raw LC-MS data files converted by MSConvert to mzXml be converted
#' as MS/MS files (e.g. data-dependent MS/MS files).
#' @param zeroFillvalue value to fill missing and zero values in the xcms peak table. If argument is left blank
#'  the default is to fill with half the smallest non-zero value. (see: \code{\link{zeroFill}}).
#' @param normMethod normalization method to perform (see: \code{\link{signNorm}}). 
#' If argument is left blank no normalization is performed. Options include 
#' median fold change (probabilistic quotient) 'medFC' or total ion signal 
#' 'totIon' methods.
#' @param manBatchAdj character vector of co-variate table column names. If argument is not supplied automatic cluster identification/ batch adjustment using the \code{\link{pcaClustId}} function will occur prior to PCA analysis. If multiple
#' column names are supplied a multiple linear regression will be calculated and the batch adjusted residuals obtained from the model.
#' @param replicates logical (default = FALSE) if TRUE then the 3rd column of the co-variate table supplied will be used to identify analytical/ preparative replicates of the same sample. This information will be used to average signal intensities of analytical replicates.
#' @param LogTransBase numeric base value for log transformation. defaults to exponential of 1 (see: \code{\link{logTrans}}).
#' @param blankFC numeric minimum fold difference for each xcms peak table feature 
#' between the sample injections (numerator) and negative control blank injections (denominator). 
#' This will result in all xcms peak table features which are less than blankFC 
#' fold change higher in the samples than the blanks being removed (default = 2 fold).
#' This background substraction will not occur until there are at least 1 blank and 1 sample
#' data file acquired and converted to mzXML files. 
#' @param pAdjustMethod character p-value multiple testing adjustment method (see: \code{\link{p.adjust}}).  
#' @param corrThresh correlation coefficient threshold (non-parametric Spearman's 
#' Rho) to group features within a retention time cluster.
#' @param minFeat minimum number of features with a Retention time/ correlation
#' cluster to consider it a group (default = 2, i.e. a cluster must contain at
#' least 2 features to be considered a group).
#' @param hclustMethod hierarchical clustering method to \code{\link{hclust.vector}} 
#' method of fastcluster package (default = "median").
#' @param distMeas distance measure for retention time clustering (default = "euclidean"). 
#' see \code{\link{hclust.vector}}.
#' Within retention clusters dissimilarity is computed as 1-correlation coefficient.
#' @param xcmsSetArgs list of arguments of the xcms function \code{\link{xcmsSet}}.
#' @param pcaOutIdArgs list of arguments to the \code{\link{pcaOutId}} function.
#' @param maxTime maximum time (in minutes) from the time the last raw data file
#' was written. If most recent raw data file is older than this time then simExTargId
#' will stop. This is designed to stop the process if necessary after an extended period
#' of time. default = 60 mins.
#' 
#' @details The function is designed to faciliate 
#' simultaneous collection of raw metabolomic
#' profiling data, conversion from an instrument manufacturer's proprietary format to the
#'  open file format mzXML conversion (using the command line
#' version of MSConvert) to a new analysis directory, xcms based peak picking/
#' alignment, data pre-processing (\code{\link{zeroFill}}, \code{\link{logTrans}},
#'  \code{\link{signNorm}}), automatic PCA-based outlier removal, scores plot
#'  cluster identification and automatic potential batch adjustment (\code{\link{pcaOutId}}, 
#'  \code{\link{pcaClustId}}, \code{\link{batchAdj}}), automatic co-variate based
#'  statistical analysis and feature deconvolution (\code{\link{coVarStatType}},
#'   \code{\link{rtCorrClust}}). The data output of the workflow can be 
#'   visualized by a shiny application (\code{\link{shiny}}) at any time 
#'   (\code{\link{targIdShiny}}) and MSn fragmentation targets rapidly identified. 
#' 
#'  The function works according to the following process:
#' 
#'  1. The function can be run before the collection of the first raw data file
#'  or after collection of any number of raw MS data files. 
#'  
#'  2. The directory location where the raw data are being written to must be provided
#'  (rawDir), as well as a comma delimited text file (.csv) containing two named columns.
#'  The first column must contain the precise names of each raw MS data file
#'  that will be eventually created in the raw data directory. 
#'  The second column of the .csv file must contain the binary class membership of
#'  each raw MS data file to be considered for statistical analysis. e.g. case/ control,
#'  treated/ untreated, exposed/ unexposed, blank/ sample 
#'  
#'  N.B. Raw data files such as blank and quality control for example which are written to 
#'  the raw data directory during an experiment will also be included in the
#'  xcms peak-picking and alignment process but will not be considered in the 
#'  subsequent statistical analysis.
#'  
#'  This information will be used to inform the 
#'  non-parametric statistical analysis at present only a binary case is supported
#'  i.e. two sample classes.    
#' 
#' 
#' analysis. 
#' 
#' @examples 
#' # example XCMS peak picking/ grouping/ alignment settings for high-resolution LC-MS data (Q-ToF). 
#' # peakwidth=c(2, 20), ppm=10, snthresh=5, bw=2, mzwid=0.015, minfrac=0.25, method="centWave"
#' @seealso \code{\link{xcmsSet}}, \code{\link{retcor}}, \code{\link{group}}, 
#' \code{\link{fillPeaks}}, \code{\link{diffreport}},
#' @export
simExTargId <- function(rawDir=NULL, analysisDir=NULL, coVar=NULL, nSlaves=NULL,
                        minFiles=10, centroid=TRUE, MS2=FALSE, 
                        zeroFillvalue=NULL, normMethod=NULL,  
                        manBatchAdj=NULL, LogTransBase=exp(1), blankFC=2, 
                        replicates=FALSE, pAdjustMethod="BH", 
                        corrThresh=0.8, minFeat=1, hclustMethod="median", 
                        distMeas="euclidean", 
                        bw=2, mzwid=0.015, minfrac=0.25,
                        xcmsSetArgs=list(peakwidth=c(2, 20), ppm=10, snthresh=5, 
                                         method="centWave"),
                        pcaOutIdArgs=list(cv="q2", scale="pareto", centre=T),
                        maxTime=60){
  
# if necessary select raw data directory  
  if(is.null(rawDir)){
    message("Select the directory your raw data files will be/ have already started to be written to...\n")
    flush.console()
    
  rawDir <- tcltk::tk_choose.dir(default = "", 
                                 caption = "Select the directory your raw data files will be written to...")
  }
# if necessary select analysis (i.e. results output) directory 
  if(is.null(analysisDir)){
    message("Select the directory your results output and mzXML files will be written to...\n")
    flush.console()
      
    analysisDir <- tcltk::tk_choose.dir(default = "", 
                                        caption = "Select the directory your .mzXML files/ results files will be written to...")
  }
  
  # create sub-directories to save mzXML files and output
  # mzXML files
  mzXmlDir <- paste0(analysisDir, "/mzXmlFiles")
  suppressWarnings(dir.create(mzXmlDir))
  # R output directory for RData files and scripts
  rDir <- paste0(analysisDir, "/R")
  suppressWarnings(dir.create(rDir))
  
  # Rdata file directory within R directory
  rDataDir <- paste0(rDir, "/rDataFiles")
  suppressWarnings(dir.create(rDataDir))
  
  # output directory for results output
  outputDir <- paste0(analysisDir, "/output")
  suppressWarnings(dir.create(outputDir))
  
  # output directory for results output
  peakTableDir <- paste0(outputDir, "/01.peakTables")
  suppressWarnings(dir.create(peakTableDir))
  # save peak picking parameters for records
  parameters.tmp <- as.data.frame(xcmsSetArgs)
  write.csv(parameters.tmp, paste0(peakTableDir, "/xcmsParameters.csv"), 
            row.names=F)
  
  # preprocessing results
  preProcDir <- paste0(outputDir, "/02.preProc")
  suppressWarnings(dir.create(preProcDir))
  # save Preprocessing parameters used for records
  parameters.tmp <- data.frame(zeroFillvalue=ifelse(is.null(zeroFillvalue), 
                                               "halfMinNonZero", zeroFillvalue),
                               logTransBase=LogTransBase, 
                               normMethod=ifelse(is.null(normMethod), 
                                                 "NoNorm", normMethod))
  write.csv(parameters.tmp, paste0(preProcDir, "/preProcParameters.csv"), row.names=F)

  # PCA results
  pcaDir <- paste0(outputDir, "/03.PCA")
  suppressWarnings(dir.create(pcaDir))
  # save Pca parameters used for records
  # extract formals
  tmpArgs <- unlist(formals(pcaOutId), recursive = T)
  # remove ellipsis
  tmpArgs <- tmpArgs[names(tmpArgs) != "..."]
  # remove any args already supplied
  argsLogic <- names(tmpArgs) %in% names(pcaOutIdArgs)  
  parameters.tmp <- as.data.frame(c(pcaOutIdArgs, tmpArgs[argsLogic == F]))
  write.csv(parameters.tmp, paste0(pcaDir, "/pcaParameters.csv"), row.names=F)
  
  # statistical analysis results
  statsDir <- paste0(outputDir, "/04.stats")
  suppressWarnings(dir.create(statsDir))
  # MS/MS targets directory
  ms2targsDir <- paste0(outputDir, "/05.MS2targs")
  suppressWarnings(dir.create(ms2targsDir))
  # documents directory for output of Rmarkdown reports and future manuscript writing
  docDir <- paste0(analysisDir, "/doc")
  suppressWarnings(dir.create(docDir))
 
  # select covariates
  if(is.null(coVar)){
    message("Select your sample runorder/ covariates file...\n")
    flush.console()
    
    coVar <- tcltk::tclvalue(tcltk::tkgetOpenFile(filetypes = "{{comma delimited text file} {.csv}} {{All files} *}",
                                                title="Select your sample runorder/ covariates file",
                                                initialdir=analysisDir))
   } 
  # read in coVar table if necessary
  coVar <- read.csv(coVar, stringsAsFactors=F, header=T)
  # extract sample file names
  sampleFileNames <- as.character(coVar[, 1])
  # extract sample classes identities
  xcmsGroups <- as.character(coVar[, 2])
  # extract potential replicates identities column
  if(ncol(coVar) > 2){
  potentialReps <- as.factor(coVar[, 3])
  }
  
  # if necessary change file names
  coVar[, 1] <- gsub('-', '.', coVar[, 1])
  # if the sample begins with a number add an X as these will appear with X 
  # in data.frame
  numBeginning <- grepl("^[0-9]", coVar[, 1])
  coVar[, 1][numBeginning] <- paste0("X", coVar[, 1][numBeginning])
  
  # detect raw samples in directory
  rawFiles <- dir(rawDir, full.names=T, pattern="\\.d$|\\.RAW$")
  
  if(length(rawFiles) == 0){
    # wait for first file to be acquired
    message("waiting for the first raw data file to be acquired...\n")
    flush.console()
    
  # while loop to detect raw files
  # once 2nd raw file created progress to peak picking
  while(length(rawFiles) == 0){
    rawFiles <- dir(rawDir, full.names=T, pattern="\\.d$|\\.RAW$")
  if(length(rawFiles) > 0){
      message("First raw data file detected...\n")
      flush.console()
      break
    }
  }
  } else {
    # wait for first file to be acquired
    message(length(rawFiles), " raw data file(s) detected...\n")
    flush.console() 
  }
  # mzXML file detect
  mzXmlFiles <- dir(mzXmlDir, full.names=T, pattern="\\.mzXML$")
  # have files converted properly (i.e.) file size greater 100 bytes
  fileSizes <- file.info(mzXmlFiles)$size > 100
  # inform user of files that may not have converted properly
  if(any(fileSizes == F)){
    message('The following mzXML files have a file size of less than 100 bytes',
            ' and may not have converted properly: \n', 
            paste0(basename(mzXmlFiles)[fileSizes == F], '\n'))
  }
  mzXmlFiles <- mzXmlFiles[fileSizes]
  # basenames of both mzXML files and raw files
  bnMzXmlFiles <- gsub("\\.mzXML$", "", basename(mzXmlFiles))
  bnRawFiles <- gsub("\\.d$|\\.RAW$", "", basename(rawFiles))
  # subset last change time at least 5 mins
  # establish last file modification time
  timeChIndx <- sapply(rawFiles, function(rawFile){
    # if the raw file is a directory then identify all files within
    dirFilesTmp <- list.files(path=rawFile, recursive=T, full.names=T)
    if(length(dirFilesTmp) > 0){
      # if any files less than 5 minutes old don't do anything as still writing file
    ready <- difftime(Sys.time(), file.info(dirFilesTmp)$mtime, units='mins') < 5
    ready <- any(ready == T)
    ready <- ifelse(ready, F, T)
    } else {
    ready <- difftime(Sys.time(), file.info(rawFile)$mtime, units='mins') > 5 
    }
    return(ready)})
  # subtract mzXmlfiles from raw
  rawFiles <- rawFiles[(bnRawFiles %in% bnMzXmlFiles) == F & timeChIndx]
  
  # compare sample file names to .mzXmL files
  SampleConv <- sampleFileNames %in% gsub("\\.mzXML$", "", basename(mzXmlFiles))
  
  # Establish msconvert commands to be sent to shell commands, centroid/ MS2 
  # (i.e. data dependent/independent) if necessary
  if(centroid == T){
    convType <- ifelse(MS2 == T, '" --32 --mzXML --filter "peakPicking true 1-2" -o "',
                       '" --32 --mzXML --filter "peakPicking true 1-" -o "')
  } else {
    convType <- ifelse(MS2 == T, '" --32 --mzXML --filter "msLevel 1-2" -o "',
                       '" --32 --mzXML --filter "msLevel 1-" -o "')
  }
  # if all already converted do nothing
  if(any(SampleConv == F)){ 
    message("Waiting for next raw data file...\n")
    flush.console()
    # if any of the samples have not been converted to mzXML then convert
  while(any(SampleConv == F)){  
    # raw files detect
    rawFiles <- dir(rawDir, full.names=T, pattern="\\.d$|\\.RAW$")
    # mzXML file detect
    mzXmlFiles <- dir(mzXmlDir, full.names=T, pattern="\\.mzXML$")
    # have files converted properly (i.e.) file size greater 100 bytes
    fileSizes <- file.info(mzXmlFiles)$size > 100
    # inform user of files that may not have converted properly
    if(any(fileSizes == F)){
      message('The following mzXML files have a file size of less than 100 bytes',
              ' and may not have converted properly: \n', 
              paste0(basename(mzXmlFiles)[fileSizes == F], '\n'))
    }
    mzXmlFiles <- mzXmlFiles[fileSizes]
    # basenames of both mzXML files and raw files
    bnMzXmlFiles <- gsub("\\.mzXML$", "", basename(mzXmlFiles))
    bnRawFiles <- gsub("\\.d$|\\.RAW$", "", basename(rawFiles))
    # subset last change time at least 5 mins
    # establish last file modification time
    timeChIndx <- sapply(rawFiles, function(rawFile){
      # if the raw file is a directory then identify all files within
      dirFilesTmp <- list.files(path=rawFile, recursive=T, full.names=T)
      if(length(dirFilesTmp) > 0){
          ready <- difftime(Sys.time(), file.info(dirFilesTmp)$mtime, units='mins')
          ready <- min(ready)
      } else {
        ready <- difftime(Sys.time(), file.info(rawFile)$mtime, units='mins')
      }
      return(ready)})
  
      # if longer than max time then stop the process.
      if(min(timeChIndx) > maxTime){
        stop("The last raw data file was written longer than an hour ago...stopping the simExTargId process")
      }
    
    # convert TimeChanged index into logical 
    # if any files less than 5 minutes old don't do anything as still writing file
    timeChIndx <- timeChIndx > 5
    # subtract mzXmlfiles from raw
    rawFiles <- rawFiles[(bnRawFiles %in% bnMzXmlFiles) == F & timeChIndx]
    # compare sample file names to .mzXmL files
    SampleConv <- sampleFileNames %in% gsub("\\.mzXML$", "", basename(mzXmlFiles))
    
    if(length(rawFiles) > 0){  
    message("Converting to mzXML and initial peak picking: \n", 
            paste(basename(rawFiles), collapse="\n"), "\n") 
    flush.console()
    
    command <- paste0('msconvert "', rawFiles, convType, mzXmlDir, '"')
    # if nSlaves not null then convert to mzXML in parallel
    if(!is.null(nSlaves) & length(command) > 1){
    # testing...
    # pmt <- proc.time()
    # start cluster
    message(paste0("Starting SNOW cluster with ", nSlaves, " local sockets...\n"))
    flush.console()
    cl <- parallel::makeCluster(nSlaves) 
    doSNOW::registerDoSNOW(cl)
    message("Converting raw files and saving in mzXmlFiles directory...\n")
    flush.console()
    # foreach and dopar from foreach package
    outTmp <- foreach(rawF=1:length(command)) %dopar% {system(command[rawF])}
    # stop SNOW cluster
    parallel::stopCluster(cl)
    rm(outTmp)
    # 6 cores 1519.29 seconds, 25.32 mins or 22.34 seconds per .RAW file 
    # (average 280 MB per .RAW file)
    # convTimeParallel <- proc.time() - pmt
    } else {    
    # single threaded
    #pmt <- proc.time()
    sapply(command, system)
    # single threaded 2500.63 seconds, 41.67 mins or 36.77 seconds per .RAW file 
    # (average 280 MB per .RAW file)
    #convTimeSingleThread <- proc.time() - pmt
    }
    message("...Done")
    flush.console()
    
    # newly converted file names
    tmpMzXmlNames <- paste0(mzXmlDir, "/", gsub("\\.d$|\\.RAW$", ".mzXML", 
                                                basename(rawFiles)))
    # identify xcmsSet.RData file in directory
    xcmsSet_file <- dir(rDataDir, full.names=T, pattern="01\\.xcmsSet\\.RData")
    # if xcmsSet R data file exists then load
    if(length(xcmsSet_file) == 1){
      load(xcmsSet_file)
      # check to see if any newly converted mzXML files are already in the
      # xcmsSet object, if so then split and remove the file
      alreadyPPindx <- rownames(tmpXcmsSet@phenoData) %in% gsub('\\.mzXML$', '', basename(tmpMzXmlNames))
      # split to remove
      tmpXcmsSet <- split(tmpXcmsSet, f=alreadyPPindx)[[1]]
      # alternatively if a file is not present in the xcmsSet object but there is a
      # mzXML file then add the file to the tmpMzXmlNames vector
     } else {
      # create empty object for conditional
      tmpXcmsSet <- list()
    }
    # perform xcms peak picking and save as RData file if no additional cores then
    # loop through mzXMLfiles
    if(is.null(nSlaves) | length(tmpMzXmlNames) == 1){
          for(mzXMLfile in tmpMzXmlNames){
                message("Detecting features in file : ", basename(mzXMLfile))
                flush.console()
                if(file.exists(mzXMLfile)){
                xcmsSetArgs$files <- mzXMLfile
                xcmsSet_tmp <- do.call(xcms::xcmsSet, xcmsSetArgs)
                # if xcmsSet R data file loaded/ tmpXcmsSet exists then concatenate
                if(length(tmpXcmsSet) > 0){
                tmpXcmsSet <- c(tmpXcmsSet, xcmsSet_tmp) 
                # if no R data file already in the directory
                } else {
                tmpXcmsSet <- xcmsSet_tmp
                } 
                } else {
                  message(mzXMLfile, '\n does not exist and may not have converted properly...')
                  flush.console()
                }
              } # end for loop   
      # else if nSlaves have been specified and more than one tmpMzXmlFile then 
      # peak pick in parallel
      } else { 
        # check to see if files have converted properly
        mzXmlFileExists <- file.exists(tmpMzXmlNames)
        if(any(mzXmlFileExists == F)){
        message(paste0(tmpMzXmlNames[mzXmlFileExists == F], "\n"), ' mzXML file(s) do not exist and may not have converted properly...')
        flush.console()
        # subset 
        tmpMzXmlNames <- tmpMzXmlNames[mzXmlFileExists]
        }
        # add temporary mzXML file names into xcmsSet arguments list  
      xcmsSetArgs$files <- tmpMzXmlNames
      xcmsSet_tmp <- do.call(xcms::xcmsSet, c(xcmsSetArgs, nSlaves=nSlaves))
      # if xcmsSet R data file loaded/ tmpXcmsSet exists then concatenate
        if(length(tmpXcmsSet) > 0){
          tmpXcmsSet <- c(tmpXcmsSet, xcmsSet_tmp) 
          # if no R data file already in the directory
        } else {
          tmpXcmsSet <- xcmsSet_tmp
        } 
      }
  # add xcmsGroup information into xcmsSet object before saving
  sampleIndx <- match(row.names(tmpXcmsSet@phenoData), sampleFileNames)
  xcms::sampclass(tmpXcmsSet) <- xcmsGroups[sampleIndx]
  # save xcmsSet class object as RData
  message("Saving 01.xcmsSet.RData file...\n")
  flush.console()
  save(tmpXcmsSet, file=paste0(rDataDir, "/01.xcmsSet.RData"))
 
  # conduct rest of xcms peak picking and stat analysis if sufficient samples 
  # acquired and picked
  if(length(sampleIndx) >= minFiles){ 
            
      # further steps of XCMS peak picking
      message("obiwarp retention time correction... \n")
      flush.console()
      
      tmpXcmsSet <- xcms::retcor(tmpXcmsSet, method="obiwarp")
      # save retcor xcmsSet class object as RData
      message("Saving 02.retcor.RData file...\n")
      flush.console()
      
      save(tmpXcmsSet, file=paste0(rDataDir, "/02.retcor.RData"))  
      
      message("Peak grouping... \n")
      flush.console()
      
      tmpXcmsSet <- xcms::group(tmpXcmsSet, bw=bw, mzwid=mzwid, minfrac=minfrac)
      # save grouped xcmsSet class object as RData
      message("Saving 03.group.RData file...\n")
      flush.console()
      
      save(tmpXcmsSet, file=paste0(rDataDir, "/03.group.RData"))  
      
      message("zero filling... \n")
      flush.console()
      
      tmpXcmsSet <- suppressWarnings(xcms::fillPeaks.chrom(tmpXcmsSet, 
                                                           nSlaves=ifelse(!is.null(nSlaves), 
                                                                          nSlaves, 0)))
      # save zerofilled xcmsSet class object as RData
      message("Saving 04.fillPeaks.RData file...\n")
      flush.console()
      
      save(tmpXcmsSet, file=paste0(rDataDir, "/04.fillPeaks.RData"))  
      
      message("01. Saving xcms peak table to 01.peakTables directory... \n")
      flush.console()
    
      # extract matrix of peak values
      xcmsPeaks <- data.frame(EICno=seq(1, nrow(tmpXcmsSet@groups), 1), 
                              peakTable(tmpXcmsSet), stringsAsFactors=F)
      colnames(xcmsPeaks)[c(2, 5)] <- c('mzmed', 'rtmed')
      # padded number of samples for file naming
      nSamplesTmp <- sprintf("%03d", length(sampleIndx))
      # write csv of current peak table if the write attempt does not work then skip and warn user
      writeCsvAttempt(df=xcmsPeaks, fileName=paste0(peakTableDir, "/xcmsOutput.csv"))

#       message("identifying near-zero-variance peak groups...\n")
#       flush.console()
      
#       nzv.tmp <- caret::nearZeroVar(xcms::groupval(tmpXcmsSet), saveMetrics=T)
#       cvAll <- apply(xcms::groupval(tmpXcmsSet), 1, function(x) (sd(x)/mean(x)) * 100)
      # deconvolute xcms output
  
    message("Preprocessing, PCA outlier removal, statistical analysis and MS2 target identification...\n\n")
    flush.console()
    # begin preprocessing of data
    message("02. Preprocessing peak table...\n")
    flush.console()
    # observation (sample) names
    obsNames <- gsub('-', '.', row.names(tmpXcmsSet@phenoData))
    # if the sample begins with a number add an X as these will appear with X 
    # in data.frame
    numBeginning <- grepl("^[0-9]", obsNames)
    obsNames[numBeginning] <- paste0("X", obsNames[numBeginning])
  
    # zero fill
    xcmsPeaks <- MetMSLine::zeroFill(xcmsPeaks, obsNames, value=zeroFillvalue)
    # subtract background i.e. blanks
    sampleGroups <- coVar[sampleIndx, 2]
    blanksIndx <- grep('blank', sampleGroups, ignore.case=T)
    samplesIndxTmp <- grep('sample', sampleGroups, ignore.case=T)
    # if there are both blanks and samples then substract the background
    if(length(blanksIndx) > 0 & length(samplesIndxTmp) > 0){
      message('Subtracting background features from peak table using blank injections: \n')
      flush.console()
      message('Any features less than ', blankFC, ' fold higher in the samples than the blanks will be removed...\n')
      flush.console()
      message(nrow(xcmsPeaks), ' xcms features were detected before background subtraction...\n')
      flush.console()
      # background subtraction fold change blanks vs. mean all samples
      xcmsPeaks$fcBlankSamples <- apply(xcmsPeaks[, obsNames[samplesIndxTmp]], 1, mean) / apply(xcmsPeaks[, obsNames[blanksIndx]], 1, mean)
      # subset peak table
      xcmsPeaks <- xcmsPeaks[xcmsPeaks$fcBlankSample > blankFC, ]
      message(nrow(xcmsPeaks), ' xcms features remain after background subtraction...\n')
      flush.console()
    }
    # log transform 
    xcmsPeaks <- MetMSLine::logTrans(xcmsPeaks, obsNames, base=LogTransBase)
    # normalization if required
    if(!is.null(normMethod)){
      xcmsPeaks <- MetMSLine::signNorm(xcmsPeaks, obsNames, method=normMethod)  
    }
    # save results
    message("saving pre-processed peak table to 02.preProc directory...\n")
    flush.console()
    # write pre-processed peak table to file
    writeCsvAttempt(df=xcmsPeaks, fileName=paste0(preProcDir, "/preProcessed.csv"))
       
    message("03. PCA projection and score cluster/ batch effect identification...\n")
    flush.console()
    # create a PCA results sub-directory
    tmpPcaResultsDir <- paste0(pcaDir, "/", nSamplesTmp, "_samples")
    suppressWarnings(dir.create(tmpPcaResultsDir))
    # if the replicates argument is TRUE then combine intensities of replicate
    # injections
    if(replicates == T){
    # if necessary change file names
    tmpObsNames <- gsub('-', '.', sampleFileNames)
    # if the sample begins with a number add an X as these will appear with X 
    # in data.frame
    numBeginning <- grepl("^[0-9]", tmpObsNames)
    tmpObsNames[numBeginning] <- paste0("X", tmpObsNames[numBeginning])
    # identify coVar and xcms peak column names
    obsIndx <- match(tmpObsNames, colnames(xcmsPeaks)) 
    obsTable <- xcmsPeaks[, obsIndx]
    # row-wise apply duplicate signal averaging
    message("replicate sample signal averaging...\n")
    flush.console()
    repAveraged.df <- t(apply(obsTable, 1, function(Var){
                        return(tapply(Var, potentialReps, mean))}))
    repAveraged.df <- as.data.frame(repAveraged.df)
    # create new obsNames
    obsNames <- tapply(sampleFileNames, potentialReps, paste, collapse="_")
    obsNames <- paste0("mean_", obsNames)
    colnames(repAveraged.df) <- obsNames
    # replace xcms peak table columns appropriately
    xcmsPeaks[, obsIndx[1:ncol(repAveraged.df)]] <- repAveraged.df
    # remove unnecessary columns 
    xcmsPeaks <- xcmsPeaks[, -obsIndx[(ncol(repAveraged.df) + 1):length(obsIndx)]]
    # add column names
    colnames(xcmsPeaks)[obsIndx[1:ncol(repAveraged.df)]] <- obsNames
    # subset coVariates table
    coVar <- coVar[duplicated(coVar[, 3]) == F, ]
    # sort table 
    coVar <- coVar[order(coVar[, 3]), ]
    # add in new names
    coVar[, 1] <- obsNames
    }
    # if manBatchAdj column names supplied then manually adjust for batch
    # prior to PCA analysis
    if(!is.null(manBatchAdj)){
    
    missingValIndx <- (is.na(coVar[, manBatchAdj]) | coVar[, manBatchAdj] == "") == F 
    coVarIndxTmp <- match(obsNames, coVar[missingValIndx, 1])
    obsIndxTmp <- which(!is.na(coVarIndxTmp))
    coVarIndxTmp <- coVar[, 1] %in% obsNames[obsIndxTmp]
    # adjust the xcmsPeaks by the manually supplied co-variate(s)
    batchAdjusted <- MetMSLine::batchAdj(xcmsPeaks, obsNames[obsIndxTmp], 
                                         coVar[coVarIndxTmp, manBatchAdj])[[2]]
    # write batchAdjusted residual xcms peak table to a .csv
    writeCsvAttempt(df=xcmsPeaks, fileName=paste0(preProcDir, "/manualBatchAdjResid_", 
                                                  paste0(manBatchAdj, collapse="_"), ".csv"))
    } else {
    batchAdjusted <- data.frame()
#   pcaOutIdArgs$obsNames <- obsNames 
    }
    # if necessary subset observation arguments
    pcaOutIdArgs$obsNames <- obsNames
    # add xcmsPeaks to pcaOutIdArgs
    pcaOutIdArgs$peakTable <- xcmsPeaks
   
    # identify potential analytical outliers
    pcaOutResults <- do.call(MetMSLine::pcaOutId, pcaOutIdArgs)
    # save results
    message("saving PCA results to 03.PCA/", basename(tmpPcaResultsDir), 
            " directory...\n")
    flush.console()
    # number of pca iterations
    nIterPca <- length(pcaOutResults$pcaResults)
    # create number for PCA indx
    xcmsGroupScoreIndx <- as.numeric(as.factor(coVar[, 2]))
    # save pca plots in results directory
    for(iter in 1:nIterPca){
    # create PNG graphics device
    png(paste0(tmpPcaResultsDir, "/pcaOutId_scores_iter", iter, ".png"), 
        width=1200, height=1200, res=150)
    # class id from sample names
    scoreSampNames <- names(pcaOutResults$pcaResults[[iter]]$possOut)
    # from coVar
    scoreSampNames <- xcmsGroupScoreIndx[match(scoreSampNames, coVar[, 1])]
    scoreSampNames <- scoreSampNames + 2
    scoreSampNames[pcaOutResults$pcaResults[[iter]]$possOut] <- 1.5
    # plot 1st two pcs using modified plotPcsEx and colour according to outlier and class
    MetMSLine::plotPcsEx(pcaOutResults$pcaResults[[iter]]$pcaResult,
                         pcaOutResults$pcaResults[[iter]]$exHotEllipse,
                         type="scores", 
                         col=scoreSampNames + 1,
                         pch=pcaOutResults$pcaResults[[iter]]$possOut + 16,
                         cex=pcaOutResults$pcaResults[[iter]]$possOut + 1)
    # close PNG graphics device
    dev.off()
    # save scores
    scoresTmp <- pcaOutResults$pcaResults[[iter]]$pcaResult@scores
    write.csv(scoresTmp, paste0(tmpPcaResultsDir, "/pcaOutId_scores_iter", 
                                iter, ".csv"))
    # save loadings
    loadingsTmp <- pcaOutResults$pcaResults[[iter]]$pcaResult@loadings
    write.csv(loadingsTmp, paste0(tmpPcaResultsDir, "/pcaOutId_loadings_iter", 
                                iter, ".csv"))
    }
    # remaining sample names
    remNames <- names(pcaOutResults$pcaResults[[nIterPca]]$possOut)
    outliers <- data.frame(outliers=setdiff(obsNames, remNames))
    # if necessary extract outlier removed peakTable
    # if xcmsPeaks and outlier removed peak table identical then do not save 
    # any new peak table
    if(nrow(outliers) != 0){
    # write outlier removed peak table to file
    xcmsPeaks <- pcaOutResults$outRem
    write.csv(outliers, paste0(tmpPcaResultsDir, "/", 
                                length(obsNames) - length(remNames), 
                                "_outliers_", nSamplesTmp, 
                                "_samples.csv"), row.names=F)
    } else {
    message("No outliers were removed therefore no outlier removed peak table will be saved in the 03.PCA/", 
            basename(tmpPcaResultsDir), " sub-directory...\n")
    flush.console()
    }
    # subset covariates table  
    remNamesCoVar <- coVar[, 1] %in% remNames
    # subset Covariates table
    coVarOutRem <- coVar[remNamesCoVar, , drop=F]
    # stats analysis
    statResults <- vector('list', ncol(coVarOutRem) - 1)
    statResultsBatchAdj <- vector('list', ncol(coVarOutRem) - 1)
    
    message("04. statistical analysis, ", length(statResults), " co-variates...\n")
    flush.console()
    # id covariates with missing or only one value
    missingIndx <- apply(coVarOutRem, 2, function(Var){ 
    any(length(unique(Var[complete.cases(Var)])) == 1)})
    if(any(missingIndx)){
      cat("The following coVarTable columns contain only one value and will not be considered:\n",
          paste0(seq(1, sum(missingIndx), 1), ". ", colnames(coVarOutRem)[missingIndx], "\n"))
      if(all(missingIndx)){
        stop("All coVarOutRem table columns contain missing values")
      }
      # subset and remove missing only one value columns
      coVarOutRem <- coVarOutRem[, missingIndx == F, drop=F]  
    }
    # create logical data frame of missing values
    coVarOutRem_logi <- data.frame(matrix(T, nrow=nrow(coVarOutRem), 
                                          ncol=ncol(coVarOutRem)))
    coVarOutRem_logi[coVarOutRem == ""] <- F
    coVarOutRem_logi[is.na(coVarOutRem)] <- F
    # if nSlaves then run in parallel
    if(!is.null(nSlaves)){
      
      # covariate name
      message("co-variate table columns:\n", 
              paste0(colnames(coVarOutRem)[2:ncol(coVarOutRem)], "\n"), "\n", 
              nrow(xcmsPeaks), " multiple comparisons\n")
      flush.console()
      
      # start cluster
      message(paste0("Starting SNOW cluster with ", nSlaves, " local sockets...\n"))
      flush.console()
      cl <- parallel::makeCluster(nSlaves) 
      doSNOW::registerDoSNOW(cl)
      
      # foreach and dopar from foreach package
      statResults <- foreach(i=2:ncol(coVarOutRem)) %dopar% {
        MetMSLine::coVarTypeStat(peakTable=xcmsPeaks, 
                                 obsNames=coVarOutRem[coVarOutRem_logi[, i], 1], 
                                 coVariate=coVarOutRem[coVarOutRem_logi[, i], i],
                                 base=LogTransBase, 
                                 MTC=pAdjustMethod)
      }
      # stop SNOW cluster
      parallel::stopCluster(cl) 
      # if batch adjustment was necessary then stats analysis on other data frame
      if(length(batchAdjusted) > 0){
        message("Now performing statistical analyses with batch adjusted data...\n")
        # start cluster
        message(paste0("Starting SNOW cluster with ", nSlaves, " local sockets...\n"))
        flush.console()
        cl <- parallel::makeCluster(nSlaves) 
        doSNOW::registerDoSNOW(cl)
        
        # foreach and dopar from foreach package
        statResultsBatchAdj <- foreach(i=2:ncol(coVarOutRem)) %dopar% {
          MetMSLine::coVarTypeStat(peakTable=batchAdjusted, 
                                   obsNames=coVarOutRem[coVarOutRem_logi[, i], 1], 
                                   coVariate=coVarOutRem[coVarOutRem_logi[, i], i],
                                   base=LogTransBase, 
                                   MTC=pAdjustMethod)
        }
        # stop SNOW cluster
        parallel::stopCluster(cl) 
      }
    # if nSlave null for loop  
    } else {
    # create progress bar
    pb <- txtProgressBar(min=0, max=length(statResults), style=3)
    
    for(i in 2:ncol(coVarOutRem)){
      # set progress bar
      Sys.sleep(0.01)
      setTxtProgressBar(pb, i - 1) 
      message("co-variate table column: ", colnames(coVarOutRem)[i], ", ", 
              nrow(xcmsPeaks), " multiple comparisons")
      flush.console()
      # auto select coVariate type
      statResults[[i - 1]] <- MetMSLine::coVarTypeStat(xcmsPeaks, 
                                                       obsNames=coVarOutRem[coVarOutRem_logi[, i], 1], 
                                                       coVariate=coVarOutRem[coVarOutRem_logi[, i], i],
                                                       base=LogTransBase, 
                                                       MTC=pAdjustMethod)
      # if necessary then stats analysis batchadjusted data
      if(length(batchAdjusted) > 0){
        message("Batch Adjusted, co-variate table column: ", colnames(coVarOutRem)[i], ", ", 
                nrow(xcmsPeaks), " multiple comparisons")
        flush.console()
        
      statResultsBatchAdj[[i - 1]] <- MetMSLine::coVarTypeStat(batchAdjusted, 
                                                               obsNames=coVarOutRem[coVarOutRem_logi[, i], 1], 
                                                               coVariate=coVarOutRem[coVarOutRem_logi[, i], i],
                                                               base=LogTransBase, 
                                                               MTC=pAdjustMethod)
      } # cond batch adjusted
    } # end single thread for loop
    } # end parallel cond
    # saving statistical analysis objects (.RData)...
    message('saving statistical analysis objects (.RData) in RDataFiles directory...')
    flush.console()
    save(statResults, file=paste0(rDataDir, "/statResults.RData"))
    # if necessary save batch adjusted statistical test results
    if(length(batchAdjusted) > 0){
      message('saving batch adjusted statistical analysis objects (.RData) in RDataFiles directory...')
      flush.console()
      save(statResultsBatchAdj, file= paste0(rDataDir, "/statResultsBatchAdj.RData"))
    }
    
    # covar names 
    coVarNames <- colnames(coVarOutRem)[2:ncol(coVarOutRem)]
    # add stat results output to peak table 
    statResults <- do.call(cbind, lapply(1:length(statResults), function(x){
          statResTmp <- statResults[[x]]
          if(!is.character(statResTmp)){
          colNamesTmp <- paste(statResTmp$method, statResTmp$MTC, coVarNames[x], 
                               sep="_")
          resDfTmp <- statResTmp$result
          colnames(resDfTmp) <- paste0(colnames(resDfTmp), "_", colNamesTmp)
          return(resDfTmp)
          }}))
    
    # insert stats analysis into peak table
    obsIndx <- match(coVarOutRem[, 1], colnames(xcmsPeaks))
    # remaining columns
    remColsIndx <- setdiff(1:ncol(xcmsPeaks), obsIndx)
    # insert stats res and write to results directory
    statResults <- data.frame(xcmsPeaks[, remColsIndx], statResults, 
                              xcmsPeaks[, obsIndx], stringsAsFactors=F)
    # if necessary extract batch adjusted test results
    if(length(batchAdjusted) > 0){
      statResultsBatchAdj <- do.call(cbind, lapply(1:length(statResultsBatchAdj),
                                                   function(x){
          statResTmp <- statResultsBatchAdj[[x]]
          if(!is.character(statResTmp)){
          colNamesTmp <- paste("batchAdj", statResTmp$method, statResTmp$MTC, coVarNames[x], 
                               sep="_")
          resDfTmp <- statResTmp$result
          colnames(resDfTmp) <- paste0(colnames(resDfTmp), "_", colNamesTmp)
          return(resDfTmp)
          }}))
      # insert stats res and write to results directory
      statResultsBatchAdj <- data.frame(xcmsPeaks[, remColsIndx], 
                                        statResultsBatchAdj, xcmsPeaks[, obsIndx],
                                        stringsAsFactors=F)
    } # cond batch adjusted
     # deconv data using RtCorrClust
     message("deconvoluting xcms peak data by retention time and correlation clustering...\n")
     flush.console()
    
     statResults <- MetMSLine::rtCorrClust(statResults, coVarOutRem[, 1], 
                                           corrThresh=corrThresh, 
                                           minFeat=minFeat, 
                                           hclustMethod=hclustMethod, 
                                           distMeas=distMeas)
    
     message("Writing statistical results output, ", nSamplesTmp, " samples...\n")
     flush.console()
     # write csv results 
     writeCsvAttempt(df=statResults[[1]], fileName=paste0(statsDir, 
                                                          "/statResults.csv"))
     writeCsvAttempt(df=statResults[[2]], fileName=paste0(statsDir, "/statResults_wMean.csv"))
    message("...done\n")
    flush.console()
   
    if(length(batchAdjusted) > 0){
      message("deconvoluting batch adjusted xcms peak data by retention time and correlation clustering...\n")
      flush.console()
      
      statResultsBatchAdj <- MetMSLine::rtCorrClust(statResultsBatchAdj, coVarOutRem[, 1], 
                                                  corrThresh=corrThresh, 
                                                  minFeat=minFeat, 
                                                  hclustMethod=hclustMethod, 
                                                  distMeas=distMeas)
      
    message("Writing batch adjusted statistical results output, ", nSamplesTmp, " samples...\n")
    flush.console()
    # write csv results 
    writeCsvAttempt(df=statResultsBatchAdj[[1]], fileName=paste0(statsDir, 
                                                         "/statResults_BatchAdj.csv"))
    writeCsvAttempt(df=statResultsBatchAdj[[2]], fileName=paste0(statsDir, "/statResults_BatchAdj_wMean.csv"))
        
    message("...done\n")
    flush.console()
    }
    # end stats test
    } else {
    message("There are not yet sufficient samples to conduct statistical analyses...\n")
    flush.console()
    message("subsequent steps of XCMS processing and data analysis will not yet occur...\n")
    flush.console()
    } # cond if any samples in covar table have been converted/ peak picked
  } # cond if any raw files have not yet been converted/ peak picked
  # if all samples mzXML converted and peak picked then break while loop
  if(all(SampleConv == T)){
  
  message("Generating final diffreport and EICs...\n")
  flush.console()
  setwd(rDataDir)
  # generate diff report
  suppressWarnings(xcms::diffreport(tmpXcmsSet, filebase="xcmsOutput_Allsamples",
                                    eicmax=10000, eicwidth=60))
  
  message("finished...\n")
  flush.console()
  break  
  }
} # end while loop
} else { 
  message("all of the samples have already been converted to mzXML files stopping")
  flush.console()
  }# end if all samples already converted
} # end function