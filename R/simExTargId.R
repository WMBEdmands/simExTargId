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
#' @param nCores numeric number of computer cores for parallel computation.
#' @param ionMode character ionization polarity must be specific (default= NULL).
#' Must be be either "negative" or "positive" or an abbreviation starting from "neg" or "pos".
#' @param metab character or data.frame. Full path to a csv file of a list of metabolites you wish to monitor (or a data.frame). See \link{peakMonitor} for more details.
#' @param minFiles the minimum number of raw data file converted to mzXML files
#' before commencing subsequent steps of XCMS processing, pre-processing, PCA and
#' statistical analysis. Default = 10.
#' @param centroid do the raw data files need to be centroided during conversion by MSConvert
#' \url{http://proteowizard.sourceforge.net/downloads.shtml}. NB. centroiding of data
#' is necessary for utilization of the "centWave" algorithm of xcms (\code{\link{findPeaks.centWave-methods}}).
#' @param mzXml logical should raw LC-MS data files converted by MSConvert to mzXml or
#' mzML file formats (default = TRUE i.e. files are converted to the mzXml format).
#' @param zeroFillValue value to fill missing and zero values in the xcms peak table. If argument is left NULL
#'  the default is to fill with half the smallest non-zero value. (see: \code{\link{zeroFill}}).
#' @param normMethod normalization method to perform (see: \code{\link{signNorm}}).
#' If argument is left blank no normalization is performed. Options include
#' median fold change (probabilistic quotient) 'medFC' or total ion signal
#' 'totIon' methods.
#' @param manBatchAdj character vector of co-variate table column names. If argument is not supplied automatic cluster identification/ batch adjustment using the \code{\link{pcaClustId}} function will occur prior to PCA analysis. If multiple
#' column names are supplied a multiple linear regression will be calculated and the batch adjusted residuals obtained from the model.
#' @param replicates logical (default = FALSE) if TRUE then the 3rd column of the co-variate table supplied will be used to identify analytical/ preparative replicates of the same sample. This information will be used to average signal intensities of analytical replicates.
#' @param LogTransBase numeric base value for log transformation. defaults to exponential of 1 (see: \code{\link{preProc}}).
#' @param smoothSpan see \link{preProc} function of MetMSLine for details (default = 0.8).
#' @param cvThresh see \link{preProc} function of MetMSLine for details (default = 30 i.e. 30\% coefficient of variation).
#' @param blankFC numeric minimum fold difference for each xcms peak table feature
#' between the sample injections (numerator) and negative control blank injections (denominator).
#' This will result in all xcms peak table features which are less than blankFC
#' fold change higher in the samples than the blanks being removed (default = 2 fold).
#' This background substraction will not occur until there are at least 1 blank and 1 sample
#' data file acquired and converted to mzXML files.
#' @param pAdjustMethod character p-value multiple testing adjustment method (see: \code{\link{p.adjust}}).
#' @param xcmsSetArgs list of arguments of the xcms function \code{\link{xcmsSet}}.
#' @param pcaOutIdArgs list of arguments to the \code{\link{pcaOutId}} function.
#' @param maxTime maximum time (in minutes) from the time the last raw data file
#' was written. If most recent raw data file is older than this time then simExTargId
#' will stop. This is designed to stop the process if necessary after an extended period
#' of time. default = 60 mins.
#' @param mailControl List of SMTP server settings see \code{\link{sendmail}} for details.
#' Example given is for google mail.
#' @param emailAddress character vector of email address(es) from and to which to send warning email
#' that run may have stopped, QCs are outlying or signal has attenuated. (if not supplied then email notifications will not be sent)
#' see \code{\link{sendmail}}.
#' @param minFeatPseudo integer the minimum number of features for a CAMERA pseudospectrum
#' to be considered (default =2 i.e. a CAMERA pseudospectrum must consist of a
#' minimum of two LC-MS peak groups). A weighted mean (weighted by the summed peak area of all samples)
#' will be calculated for each pseudospectrum group. This removal of artefactual
#' LC-MS features is performed to reduce the multiple testing burden prior to
#' statistical analysis.
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
#'   (\code{\link{targetId}}) and MS2 fragmentation targets rapidly identified.
#'
#'  The function works according to the following process:
#'
#'  1. The function can be run before the collection of the first raw data file
#'  or after collection of any number of raw MS data files.
#'
#'  2. The directory location where the raw data are being written to must be provided
#'  (rawDir), as well as a comma delimited text file (.csv) containing at minimum
#'  3 columns namely.
#'  \enumerate{
#'  \item The first column must contain the precise names of each raw MS data file
#'  that will be eventually created in the raw data directory.
#'
#'   \item The second column of the .csv file must contain the sample class of
#'  each raw MS data file. There are 3 suggested injection/run type names that should be included
#'  in this column for Good Laboratory Practice (GLP) and ideal
#'  metabolomic experimental design. The only mandatory injection/run type name
#'  is "sample" (not case-sensitive), simExTargId will stop if this is not found in the second column.
#'
#'  This will direct simExTargId to perform statistical analysis on samples only for which
#'  co-variates are available in columns 3 -  total number of columns.
#'   In order to use the full functionality of simExTargId it is also suggested that
#'  "local" equivolume pooled quality control samples are included in the 2nd column
#'  and must be named as "QC" (not case-sensitive). If these injection types are detected then many additional
#'  processing functions and monitoring methods will be accessed. These QC
#'  samples can be monitored for signal attenuation, used to smooth the data (see
#'  \link{loessSmooth}), filter based on minimum CV\% (see \link{cvCalc}) and
#'  used to identify if the last QC injected is outlying in the automatic principle components analysis  (see \link{pcaOutId}).
#'  In the case of signal attenuation (see \link{peakMonitor}) and a QC sample
#'  being detected as an outlier an email will be sent to the email address(es)
#'  supplied.
#'
#'  In order to distinguish a column conditioning QC from a regular QC sample
#'  the prefix "cc" should be added to the name "QC". If a QC is not used for
#'  column conditioning then just "cc" (not case-sensitive) is sufficient. These
#'  will not be monitored by the peakMonitor of pcaOutId functions.
#'
#'  Finally a blank sample should be denoted as "blank" (also not case-sensitive).
#'  If these are added to the coVariate table second column then blank subtraction
#'  will also be performed (see \link{blankSub}).
#'  Any additional injection/run types are allowed e.g. "gQC" (global QC) or "MS2" for example.
#'  If the file is denoted as "MS2" then these files will be converted differently
#'  to normal MS1 profiling data by MSConvert.
#'
#'  N.B. Raw data files such as blank and quality control for example which are written to
#'  the raw data directory during an experiment will also be included in the
#'  xcms peak-picking and alignment process but will of course not be considered in the
#'  subsequent statistical analysis.
#'
#'  \item columns 3 + all subsequent columns can contain a mixture of co-variates
#'  these can be any combination of continuous or categorical variable associated
#'  with the "sample" injections/runs specified in column 2. The function
#'  \link{coVarTypeStat} will select an appropriate univariate statistical method type
#'  to use based on this contents of these columns at any stage of the profiling
#'  run. For example if 3 categorical classes of minimum size each are identified
#'  then an ANOVA analysis will be performed or if a variable is found to be
#'  continuous then a correlation analysis will be performed.
#' }
#' @examples
#' # example XCMS peak picking/ grouping/ alignment settings for high-resolution LC-MS data (Q-ToF).
#' # peakwidth=c(2, 20), ppm=10, snthresh=5, bw=2, mzwid=0.015, minfrac=0.25, method="centWave"
#' @seealso \code{\link{xcmsSet}}, \code{\link{retcor}}, \code{\link{group}},
#' \code{\link{fillPeaks}}, \code{\link{diffreport}},
#' @export
simExTargId <- function(rawDir=NULL, studyName='exampleStudyName',
                        analysisDir=NULL, coVar=NULL,
                        nCores=NULL, ionMode=NULL, metab=NULL,
                        minFiles=10, centroid=TRUE, mzXml=TRUE,
                        zeroFillValue=NULL, normMethod=NULL,
                        manBatchAdj=NULL, LogTransBase=exp(1),
                        smoothSpan=0.8, cvThresh=30, blankFC=2,
                        replicates=FALSE, pAdjustMethod="BH",

                         xcmsSetArgs=list(peakwidth=c(5, 20), ppm=10, snthresh=5,
                                         method="centWave"),
                        bw=2, mzwid=0.015, minfrac=0.25,
                        pcaOutIdArgs=list(cv="q2", scale="pareto", centre=TRUE),
                        peakMonitorArgs=list(ppm=10, rtdev=30, maxSignalAtt=20,
                                             percBelow=20),
                        maxTime=60, emailAddress=NULL,
                        mailControl=list(smtpServer="ASPMX.L.GOOGLE.COM"),
                        minFeatPseudo=2){

  if(is.null(ionMode)){
    stop('ionMode must be specified must be negative or positive or an abbreviation starting from neg or pos')
  }
  if(!grepl('neg|pos', ionMode, ignore.case = TRUE)){
    stop('ionMode argument must be be either negative or positive or an abbreviation starting from neg or pos')
  }
  if(!require(pcaMethods)){
    stop('package pcaMethods is required please install from Bioconductor')
  }

  # file conversion type (i.e. data dependent/independent) if necessary
  convTypeMS1 <- paste0('" --32 --zlib --', ifelse(mzXml, 'mzXML', 'mzML'),
                        ' --filter "', ifelse(centroid, 'peakPicking true 1-"',
                                              'msLevel 1-"'), ' -o "')

  convTypeMS2 <- paste0('" --32 --zlib --', ifelse(mzXml, 'mzXML', 'mzML'),
                        ' --filter "', ifelse(centroid, 'peakPicking true 1-2"',
                                              'msLevel 1-2"'), ' -o "')

  # email record
  emailRec <- as.character()
  # read camera rule set from ext data
  # CAMERA
  # taken from Beyond Profiling Stanstrup et. al.
  # Select which rules to use
  if(grepl('pos', ionMode, ignore.case = TRUE)){
    cameraRules <- system.file("extdata", "rules_jan_pos.csv",
                               package = "simExTargId")
  } else {
    cameraRules <- system.file("extdata", "rules_jan_neg.csv",
                               package = "simExTargId")
  }
  cameraRules <- read.csv(cameraRules, header=TRUE, stringsAsFactors = FALSE)

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

  if(is.null(metab)){
    message("Select the any metabolites you wish to monitor in the xcms output table...\n",
            'N.B. the table (.csv) must a minimum of the following mandatory columns:\n1. "name"\n2. "mzmed"\n3. "rtmed"\n')
    flush.console()
    metab <- tcltk::tclvalue(tcltk::tkgetOpenFile(title="select your metabolites to monitor file. Click cancel if you do not wish to monitor any metabolites.", initialdir=analysisDir))
  }

  if(is.character(metab)){
    if(metab != ''){
    if(!grepl('\\.csv$', metab)){
      stop('metabolite table file must be a .csv')
    }
    metab <- as.data.frame(data.table::fread(metab, header=TRUE, stringsAsFactors = FALSE))
    reqCols <- c('name', 'mzmed', 'rtmed') %in% names(metab)
    if(any(reqCols == FALSE)){
      stop('The metabolite table must consist of a minimum of 3 column names:\n1. "name"\n2. "mzmed"\n3. "rtmed"\n')
    }
    # all unique compound names
    if(any(duplicated(metab$name))){
      stop('metabolite table names must be unique. If you have the same compound name twice (i.e. with different retention times place an additional unique integer at the end of the name. e.g. lysoPC (14:0) 1, lysoPC (14:0) 1')
    }
    }
  }
  # create sub-directories to save mzXML files and output
  # mzXML files
  mzXmlDir <- paste0(analysisDir, '/', studyName, ifelse(mzXml, "_mzXmlFiles", '_mzMlFiles'))
  suppressWarnings(dir.create(mzXmlDir))
  # add mode
  mzXmlDir <- paste0(mzXmlDir, '/', toupper(substring(ionMode, 1, 3)))
  suppressWarnings(dir.create(mzXmlDir))
  # MS1 data
  mzXmlDir <- paste0(mzXmlDir, '/MS1/')
  suppressWarnings(dir.create(mzXmlDir))

  # study dir
  analysisDir <- paste0(analysisDir, '/', studyName, '_analysis')
  suppressWarnings(dir.create(analysisDir))
  analysisDir <- paste0(analysisDir, '/', toupper(substring(ionMode, 1, 3)))
  suppressWarnings(dir.create(analysisDir))
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
  parameters.tmp <- data.frame(zeroFillValue=ifelse(is.null(zeroFillValue),
                                               "halfMinNonZero", zeroFillValue),
                               logTransBase=LogTransBase,
                               normMethod=ifelse(is.null(normMethod),
                                                 "NoNorm", normMethod))
  write.csv(parameters.tmp, paste0(preProcDir, "/preProcParameters.csv"), row.names=F)

  # PCA results
  pcaDir <- paste0(outputDir, "/03.PCA")
  suppressWarnings(dir.create(pcaDir))
  # save Pca parameters used for records
  # extract formals
  tmpArgs <- unlist(formals(MetMSLine::pcaOutId), recursive = TRUE)
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
  # data directory for files that aren't modified
  dataDir <- paste0(analysisDir, "/data")
  suppressWarnings(dir.create(dataDir))

  # select covariates
  if(is.null(coVar)){
    message("Select your sample runorder/ covariates file...\n")
    flush.console()

    coVar <- tcltk::tclvalue(tcltk::tkgetOpenFile(filetypes = "{{comma delimited text file} {.csv}} {{All files} *}",
                                                title="Select your sample runorder/ covariates file",
                                                initialdir=analysisDir))
  }
  coVarBaseName <- basename(coVar)
  # read in coVar table if necessary
  coVar <- read.csv(coVar, stringsAsFactors=FALSE, header=TRUE)
  if(ncol(coVar) < 3){
    stop('worklist/co-variates table must contain a minimum of 3 columns. See system.file("extdata", "exampleWorklist_coVariates.csv", package = "simExTargId") for an example worklist/co-variates table.')
  }

  message('Saving co-variates/worklist table (.csv) to ', paste0(dataDir, '/', coVarBaseName), '\n\n')
  # write co-variates worklist to data directory
  write.csv(coVar, paste0(dataDir, '/', coVarBaseName), row.names = FALSE)
  # 1. check no missing values
  if(any(coVar[, 1] == '')){
    stop('missing file names in first column of worklist/co-variates table.  See system.file("extdata", "exampleWorklist_coVariates.csv", package = "simExTargId") for an example worklist/co-variates table.')
  }
  # 2. check all file names are unique
  dupFileNames <- duplicated(coVar[, 1])
  if(any(dupFileNames != FALSE)){
    stop('worklist/co-variates table contains non-unique file names in the first column:\n',
         paste0(coVar[dupFileNames, 1], collapse='\n'), '\n See system.file("extdata", "exampleWorklist_coVariates.csv", package = "simExTargId") for an example worklist/co-variates table.')
  }
  # 3. check for forbidden characters in file names
  forbidSym <- "/|\\\\|\\^|\\$|\\.|\\?|\\*|\\||\\+|\\(|\\)|\\[|\\{|\\-"
  if(any(grepl(forbidSym, coVar[, 1]))){
    stop('worklist/co-variates table contains forbidden symbols in file names in the first column:\n i.e. "\\\\/^$.?*|+()[{-"\n')
  }

  # 4. check sample columns are identified:
  sampTypeSumm <- table(coVar[, 2])
  if(all(grepl('sample', names(sampTypeSumm), ignore.case = TRUE) == FALSE)){
    stop('Mandatory run/injection type "sample" (not case-sensitive) was not found in the second column of the worklist/co-variates table.  See system.file("extdata", "exampleWorklist_coVariates.csv", package = "simExTargId") for an example worklist/co-variates table.')
  }

  # 5. print summary
  message('The following injection/run types were detected in column 2 of the worklist/co-variates table:\n', paste0(names(sampTypeSumm), ' (n=', sampTypeSumm, ')\n\n'), 'Is this correct?')
  flush.console()
  # proceed <- readline(paste0('\n[y/n]:'))
  #
  # if(tolower(proceed) != 'y'){
  #   warning('simExTargId run stopped. See system.file("extdata", "exampleWorklist_coVariates.csv", package = "simExTargId") for an example worklist/co-variates table.\n', immediate. = TRUE)
  # }
  # quality control idx
  qcIdx <- grepl('^QC$', coVar[, 2], ignore.case = TRUE)
  # blank idx
  blIdx <- grepl('^blank$', coVar[, 2], ignore.case = TRUE)
  # column conditioning
  ccIdx <- grepl('^cc$|^ccQC$', coVar[, 2], ignore.case = TRUE)
  # ms/ms idx
  ms2Idx <- grepl('^MS2$', coVar[, 2], ignore.case = TRUE)

  # create seperate ms2 directory
  # MS2 data
  ms2Dir <- paste0(dirname(mzXmlDir), '/MS2/')
  suppressWarnings(dir.create(ms2Dir))
  # extract sample file names
  sampleFileNames <- as.character(coVar[, 1])
  # extract sample classes identities
  sampleClasses <- as.character(coVar[, 2])
  names(sampleClasses) <- sampleFileNames
  # characterize sample classes
  sampIdx <- grepl('^sample$', sampleClasses, ignore.case = TRUE)
  # extract potential replicates identities column
  if(ncol(coVar) > 2){
  potentialReps <- as.factor(coVar[, 3])
  }

  # if necessary change file names
  #coVar[, 1] <- gsub('-', '.', coVar[, 1])
  # if the sample begins with a number add an X as these will appear with X
  # in data.frame
  numBeginning <- grepl("^[0-9]", coVar[, 1])
  if(any(numBeginning)){
    stop('File names must not begin with a number.')
  }
  # coVar[, 1][numBeginning] <- paste0("X", coVar[, 1][numBeginning])

  # detect raw samples in directory
  rawFiles <- dir(rawDir, full.names=TRUE, pattern="\\.d$|\\.RAW$")

  if(length(rawFiles) == 0){
    # wait for first file to be acquired
    message("waiting for the first raw data file to be acquired...\n")
    flush.console()

  # while loop to detect raw files
  # once 2nd raw file created progress to peak picking
  while(length(rawFiles) == 0){
    rawFiles <- dir(rawDir, full.names=TRUE, pattern="\\.d$|\\.RAW$")
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
  mzXmlFiles <- dir(mzXmlDir, full.names=TRUE, pattern="\\.mzXML$|\\.mzML$")
  # add ms2 files
  mzXmlFiles <- c(mzXmlFiles, dir(ms2Dir, full.names=TRUE, pattern="\\.mzXML$|\\.mzML$"))
  # have files converted properly (i.e.) file size greater 100 bytes
  # fileSizes <- file.info(mzXmlFiles)$size > 100
  # # inform user of files that may not have converted properly
  # if(any(fileSizes == FALSE)){
  #   message('The following mzXML/mzML files have a file size of less than 100 bytes',
  #           ' and may not have converted properly: \n',
  #           paste0(basename(mzXmlFiles)[fileSizes == F], '\n'))
  # }
  # mzXmlFiles <- mzXmlFiles[fileSizes]
  # basenames of both mzXML files and raw files
  bnMzXmlFiles <- gsub("\\.mzXML$|\\.mzML$", "", basename(mzXmlFiles))
  bnRawFiles <- gsub("\\.d$|\\.RAW$", "", basename(rawFiles))
  # subset last change time at least 5 mins
  # establish last file modification time
  timeChIndx <- sapply(rawFiles, function(rawFile){
    # if the raw file is a directory then identify all files within
    dirFilesTmp <- list.files(path=rawFile, recursive=TRUE, full.names=TRUE)
    if(length(dirFilesTmp) > 0){
      # if any files less than 5 minutes old don't do anything as still writing file
    ready <- difftime(Sys.time(), file.info(dirFilesTmp)$mtime, units='mins') < 5
    ready <- any(ready == TRUE)
    ready <- ifelse(ready, FALSE, TRUE)
    } else {
    ready <- difftime(Sys.time(), file.info(rawFile)$mtime, units='mins') > 5
    }
    return(ready)})
  # subtract mzXmlfiles from raw
  rawFiles <- rawFiles[(bnRawFiles %in% bnMzXmlFiles) == FALSE & timeChIndx]

  # compare sample file names to .mzXmL files
  SampleConv <- sampleFileNames %in% gsub("\\.mzXML$|\\.mzML$", "", basename(mzXmlFiles))

  # if all already converted do nothing
  if(any(SampleConv == FALSE)){
    message("Waiting for next raw data file...\n")
    flush.console()
    # if any of the samples have not been converted to mzXML then convert
  while(any(SampleConv == FALSE)){
    # raw files detect
    rawFiles <- dir(rawDir, full.names=TRUE, pattern="\\.d$|\\.RAW$")
    # mzXML file detect
    mzXmlFiles <- dir(mzXmlDir, full.names=TRUE, pattern="\\.mzXML$|\\.mzML$")
    # add ms2 files
    mzXmlFiles <- c(mzXmlFiles, dir(ms2Dir, full.names=TRUE, pattern="\\.mzXML$|\\.mzML$"))
    # have files converted properly (i.e.) file size greater 100 bytes
    # fileSizes <- file.info(mzXmlFiles)$size > 100
    # # inform user of files that may not have converted properly
    # if(any(fileSizes == FALSE)){
    #   message('The following mzXML files have a file size of less than 100 bytes',
    #           ' and may not have converted properly: \n',
    #           paste0(basename(mzXmlFiles)[fileSizes == FALSE], '\n'))
    # }
    # mzXmlFiles <- mzXmlFiles[fileSizes]
    # basenames of both mzXML files and raw files
    bnMzXmlFiles <- gsub("\\.mzXML$|\\.mzML$", "", basename(mzXmlFiles))
    bnRawFiles <- gsub("\\.d$|\\.RAW$", "", basename(rawFiles))
    # subset last change time at least 5 mins
    # establish last file modification time
    timeChIndx <- sapply(rawFiles, function(rawFile){
      # if the raw file is a directory then identify all files within
      dirFilesTmp <- list.files(path=rawFile, recursive=TRUE, full.names=TRUE)
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
    rawFiles <- rawFiles[(bnRawFiles %in% bnMzXmlFiles) == FALSE & timeChIndx]
    sampClTmp <- bnRawFiles[(bnRawFiles %in% bnMzXmlFiles) == FALSE & timeChIndx]
    # compare sample file names to .mzXmL files
    SampleConv <- sampleFileNames %in% gsub("\\.mzXML$|\\.mzML$", "", basename(mzXmlFiles))

    if(length(rawFiles) > 0){
    message("Converting to mzXML/mzML and initial peak picking: \n",
            paste(basename(rawFiles), collapse="\n"), "\n")
    flush.console()
    # if any are ms2 convert differently
    sampClTmp <- sampleClasses[sampClTmp]
    ms2Conv <- grepl('MS2', sampClTmp)
    command <- paste0('msconvert "', rawFiles, ifelse(ms2Conv, convTypeMS2, convTypeMS1),
                      ifelse(ms2Conv, ms2Dir, mzXmlDir), '"')
    # if nCores not null then convert to mzXML in parallel
    if(!is.null(nCores) & length(command) > 1){
    # testing...
    # pmt <- proc.time()
    # start cluster
    require(foreach)
    require(doSNOW)
    message(paste0("Starting SNOW cluster with ", nCores, " local sockets...\n"))
    flush.console()
    cl <- parallel::makeCluster(nCores, outfile='')
    doSNOW::registerDoSNOW(cl)
    message("Converting raw files and saving in mzXmlFiles directory...\n")
    flush.console()
    progress <- function(n) cat(paste0(n, ' of ', length(command),
                                       ' complete (', basename(rawFiles)[n],
                                       ').\n'))
    opts <- list(progress=progress)

    # foreach and dopar from foreach package
    outTmp <- foreach(rawF=1:length(command), .options.snow=opts) %dopar% {system(command[rawF])}
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
    tmpMzXmlNames <-  dir(mzXmlDir, full.names=TRUE, pattern="\\.mzXML$|\\.mzML$")
    # have files converted properly (i.e.) file size greater 100 bytes
    fileSizes <- file.info(tmpMzXmlNames)$size > 100
    # inform user of files that may not have converted properly
    if(any(fileSizes == FALSE)){
      message('The following mzXML files have a file size of less than 100 bytes',
              ' and may not have converted properly: \n',
              paste0(basename(tmpMzXmlNames)[fileSizes == FALSE], '\n'))
    }
    tmpMzXmlNames <- tmpMzXmlNames[fileSizes]
    # identify xcmsSet.RData file in directory
    xcmsSet_file <- dir(rDataDir, full.names=TRUE, pattern="01\\.xcmsSet\\.RData")
    # if xcmsSet R data file exists then load
    if(length(xcmsSet_file) == 1){
      load(xcmsSet_file)
      # check to see if any newly converted mzXML files are already in the
      # xcmsSet object, if so then split and remove the file
      alreadyPPindx <- rownames(tmpXcmsSet@phenoData) %in% gsub('\\.mzXML$|\\.mzML$', '', basename(tmpMzXmlNames))
      # split to remove
      tmpXcmsSet <- split(tmpXcmsSet, f=alreadyPPindx)[[1]]

      alreadyPPindx <- match(gsub('\\.mzXML$|\\.mzML$', '', basename(tmpMzXmlNames)), rownames(tmpXcmsSet@phenoData))
      tmpMzXmlNames <- tmpMzXmlNames[is.na(alreadyPPindx)]
      # alternatively if a file is not present in the xcmsSet object but there is a
      # mzXML file then add the file to the tmpMzXmlNames vector
     } else {
      # create empty object for conditional
      tmpXcmsSet <- list()
    }
    # perform xcms peak picking and save as RData file if no additional cores then
    # loop through mzXMLfiles
    if(is.null(nCores) | length(tmpMzXmlNames) == 1){
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
      # else if nCores have been specified and more than one tmpMzXmlFile then
      # peak pick in parallel
      } else {
        # check to see if files have converted properly
        mzXmlFileExists <- file.exists(tmpMzXmlNames)
        if(any(mzXmlFileExists == FALSE)){
        message(paste0(tmpMzXmlNames[mzXmlFileExists == FALSE], "\n"), ' mzXML file(s) do not exist and may not have converted properly...')
        flush.console()
        # subset
        tmpMzXmlNames <- tmpMzXmlNames[mzXmlFileExists]
        }
        # add temporary mzXML file names into xcmsSet arguments list
      xcmsSetArgs$files <- tmpMzXmlNames
      message("Detecting features in ", length(tmpMzXmlNames), " file(s)")
      flush.console()
      xcmsSet_tmp <- do.call(xcms::xcmsSet, c(xcmsSetArgs, nSlaves=nCores))
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
  sampClassTmp <- sampleClasses[sampleIndx]
  sampClassTmp <- ifelse(is.na(sampClassTmp), 'unknown', sampClassTmp)
  xcms::sampclass(tmpXcmsSet) <- sampClassTmp
  # save xcmsSet class object as RData
  message("Saving 01.xcmsSet.RData file...\n")
  flush.console()
  save(tmpXcmsSet, file=paste0(rDataDir, "/01.xcmsSet.RData"))

  # conduct rest of xcms peak picking and stat analysis if sufficient samples
  # acquired and picked
  if(sum(grepl('^sample$', sampClassTmp)) >= minFiles){

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
                                                           nSlaves=ifelse(!is.null(nCores),
                                                                          nCores, 0)))
      # save zerofilled xcmsSet class object as RData
      message("Saving 04.fillPeaks.RData file...\n")
      flush.console()

      save(tmpXcmsSet, file=paste0(rDataDir, "/04.fillPeaks.RData"))

      #Create an xsAnnotate object
      tmpXcmsSet <- CAMERA::xsAnnotate(tmpXcmsSet, polarity = ifelse(grepl('neg', ionMode, ignore.case = TRUE), 'negative', 'positive'))

      #Group after RT value of the xcms grouped peak
      tmpXcmsSet <- CAMERA::groupFWHM(tmpXcmsSet, perfwhm=1)

      #Verify grouping
      #calcCiS just because then we don't the raw data
      tmpXcmsSet <- CAMERA::groupCorr(tmpXcmsSet, graphMethod='lpc', pval=0.000001, calcCaS = TRUE,
                        cor_exp_th=0.7, cor_eic_th=0.7,calcCiS = FALSE, calcIso = FALSE)

      #Annotate isotopes
      tmpXcmsSet <- CAMERA::findIsotopes(tmpXcmsSet, ppm=10, mzabs=0.01, intval="into")

      #Annotate adducts and neutral losses
      tmpXcmsSet <- CAMERA::findAdducts(tmpXcmsSet, polarity=ifelse(grepl('neg', ionMode, ignore.case = TRUE), 'negative', 'positive'),
                           ppm=10, mzabs=0.01,
                           multiplier=4, rules=cameraRules)


      #Generate result
      message("05.CAMERA.RData file... \n")
      flush.console()

      save(tmpXcmsSet, file=paste0(rDataDir, "/05.CAMERA.RData"))

      #Generate result
      message("01. Saving xcms peak table to 01.peakTables directory... \n")
      flush.console()

      # extract matrix of peak values
      xcmsPeaks <- data.frame(CAMERA::getPeaklist(tmpXcmsSet), stringsAsFactors=FALSE)
      xcmsPeaks <- cbind(EICno=as.integer(seq(1, nrow(xcmsPeaks), 1)), xcmsPeaks)
      colnames(xcmsPeaks)[c(2, 5)] <- c('mzmed', 'rtmed')
      # make sure peak table is in run order according to coVars table
      # # observation (sample) names
      obsNames <- row.names(tmpXcmsSet@xcmsSet@phenoData)
      # sort peak table in to acquisition order if neccessary
      varCols <- setdiff(colnames(xcmsPeaks), obsNames)
      acqOrder <- match(obsNames, coVar[, 1])
      # if missing
      acqOrder <- zoo::na.spline(acqOrder)
      obsNames <- obsNames[order(acqOrder)]
      xcmsPeaks <- xcmsPeaks[, c(varCols, obsNames)]
      # save raw xcms peaks info for peak monitor after outlier removal and pre-processing
      xcmsPeaksRaw <- xcmsPeaks
      # padded number of samples for file naming
      nSamplesTmp <- sprintf("%03d", length(sampleIndx))
      # write csv of current peak table if the write attempt does not work then skip and warn user
      writeCsvAttempt(df=xcmsPeaks, fileName=paste0(peakTableDir, "/xcmsOutput.csv"))

#       message("identifying near-zero-variance peak groups...\n")
#       flush.console()

#       nzv.tmp <- caret::nearZeroVar(xcms::groupval(tmpXcmsSet), saveMetrics=TRUE)
#       cvAll <- apply(xcms::groupval(tmpXcmsSet), 1, function(x) (sd(x)/mean(x)) * 100)
      # deconvolute xcms output

    message("Preprocessing, PCA outlier removal, statistical analysis and MS2 target identification...\n\n")
    flush.console()
    # begin preprocessing of data
    message("02. Preprocessing peak table...\n")
    flush.console()

    # if the sample begins with a number add an X as these will appear with X
    # in data.frame
    # numBeginning <- grepl("^[0-9]", obsNames)
    # obsNames[numBeginning] <- paste0("X", obsNames[numBeginning])
    sampleGroups <- coVar[sampleIndx, 2]
    blanksIndxTmp <- grep('^blank$', sampleGroups, ignore.case=TRUE)
    samplesIndxTmp <- grep('^sample$', sampleGroups, ignore.case=TRUE)
    qcIndxTmp <- grep('^QC$', sampleGroups, ignore.case=TRUE)
    # minimum of 4 qcs
    qcNames <- NULL
    if(length(qcIndxTmp) >= 3){
    qcNames <- coVar[sampleIndx[qcIndxTmp], 1]
    samplesIndxTmp <- samplesIndxTmp[samplesIndxTmp < max(qcIndxTmp)]
    blanksIndxTmp <- blanksIndxTmp[blanksIndxTmp < max(qcIndxTmp)]
    }

    blankNames <- NULL
    if(length(blanksIndxTmp) > 0){
    blankNames <- coVar[sampleIndx[blanksIndxTmp], 1]
    }
    sampNames <- coVar[sampleIndx[samplesIndxTmp], 1]

    xcmsPeaks <- MetMSLine::preProc(xcmsPeaks, obsNames, sampNames,
                                    qcNames, blankNames, zeroFillValue,
                                    normMethod, cvThresh, nCores,
                                    outputDir = NULL, smoothSpan,
                                    folds = 7, baseLogT = LogTransBase,
                                    blankFCmethod = "mean",
                                    blankFCthresh = blankFC)
    # save results
    message("saving pre-processed peak table to 02.preProc directory...\n")
    flush.console()
    # write pre-processed peak table to file
    writeCsvAttempt(df=xcmsPeaks, fileName=paste0(preProcDir, "/preProcessed.csv"))

    message("03. PCA projection and score cluster/ batch effect identification...\n")
    flush.console()
    # create a PCA results sub-directory
    tmpPcaResultsDir <- paste0(pcaDir, "/", nSamplesTmp, "_files")
    suppressWarnings(dir.create(tmpPcaResultsDir))
    # if the replicates argument is TRUE then combine intensities of replicate
    # injections
    # if(replicates == TRUE){
    # # if necessary change file names
    # tmpObsNames <- gsub('-', '.', sampleFileNames)
    # # if the sample begins with a number add an X as these will appear with X
    # # in data.frame
    # # numBeginning <- grepl("^[0-9]", tmpObsNames)
    # # tmpObsNames[numBeginning] <- paste0("X", tmpObsNames[numBeginning])
    # # identify coVar and xcms peak column names
    # obsIndx <- match(tmpObsNames, colnames(xcmsPeaks))
    # obsTable <- xcmsPeaks[, obsIndx]
    # # row-wise apply duplicate signal averaging
    # message("replicate sample signal averaging...\n")
    # flush.console()
    # repAveraged.df <- t(apply(obsTable, 1, function(Var){
    #                     return(tapply(Var, potentialReps, mean))}))
    # repAveraged.df <- as.data.frame(repAveraged.df)
    # # create new obsNames
    # obsNames <- tapply(sampleFileNames, potentialReps, paste, collapse="_")
    # obsNames <- paste0("mean_", obsNames)
    # colnames(repAveraged.df) <- obsNames
    # # replace xcms peak table columns appropriately
    # xcmsPeaks[, obsIndx[1:ncol(repAveraged.df)]] <- repAveraged.df
    # # remove unnecessary columns
    # xcmsPeaks <- xcmsPeaks[, -obsIndx[(ncol(repAveraged.df) + 1):length(obsIndx)]]
    # # add column names
    # colnames(xcmsPeaks)[obsIndx[1:ncol(repAveraged.df)]] <- obsNames
    # # subset coVariates table
    # coVar <- coVar[duplicated(coVar[, 3]) == F, ]
    # # sort table
    # coVar <- coVar[order(coVar[, 3]), ]
    # # add in new names
    # coVar[, 1] <- obsNames
    # }
    # if manBatchAdj column names supplied then manually adjust for batch
    # prior to PCA analysis
    if(!is.null(manBatchAdj)){

    missingValIndx <- (is.na(coVar[, manBatchAdj]) | coVar[, manBatchAdj] == "") == FALSE
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
    # http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    # save pca plots in results directory
    for(iter in 1:nIterPca){
    # create PNG graphics device
    png(paste0(tmpPcaResultsDir, "/pcaOutId_scores_iter", iter, ".png"),
        width=1200, height=1200, res=150)
    # class id from sample names
    scoreSampNames <- names(pcaOutResults$pcaResults[[iter]]$possOut)
    # from coVar
    idxTmp <- match(scoreSampNames, coVar[, 1])
    colsPcaIt <- rep('darkgrey', length(idxTmp))
    names(colsPcaIt)[is.na(idxTmp)] <- 'unknown'
    names(colsPcaIt)[!is.na(idxTmp)] <- coVar[idxTmp[!is.na(idxTmp)], 2]
    colsPcaIt[grepl('^sample$', names(colsPcaIt))] <- 'black'
    colsPcaIt[grepl('^blank$', names(colsPcaIt))] <- 'blue'
    colsPcaIt[grepl('^QC$', names(colsPcaIt))] <- 'red'
    colsPcaIt[grepl('^ccQC$|^cc$', names(colsPcaIt))] <- 'tomato'
    remGroups <- unique(names(colsPcaIt[colsPcaIt == 'darkgrey']))
    remGroups <- remGroups[remGroups != 'unknown']
    if(length(remGroups) > 0){
      remColIdx <- match(names(colsPcaIt), remGroups)
      colsPcaIt[!is.na(remColIdx)] <- gg_color_hue(4 + length(remGroups))[remColIdx[!is.na(remColIdx)] + 4]
    }
    scoreSampNames[pcaOutResults$pcaResults[[iter]]$possOut] <- 1.5
    # plot 1st two pcs using modified plotPcsEx and colour according to outlier and class
    MetMSLine::plotPcsEx(pcaOutResults$pcaResults[[iter]]$pcaResult,
                         pcaOutResults$pcaResults[[iter]]$exHotEllipse,
                         type="scores",
                         col=colsPcaIt, #scoreSampNames + 1,
                         pch=20, #pcaOutResults$pcaResults[[iter]]$possOut + 16,
                         cex=pcaOutResults$pcaResults[[iter]]$possOut + 1,
                         main=paste0("PCA scores - ", nSamplesTmp, ' files (iteration ', iter, ')'))
    legend("topright", unique(names(colsPcaIt)), col=unique(colsPcaIt), pch = 20,
           ncol=2, text.col=unique(colsPcaIt),inset = .15)

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
      # if any outliers QC then email
      outIdx <- coVar[, 1] %in% outliers$outliers
      qcOutIdx <- grepl('^QC$', coVar[outIdx, 2])
      if(any(qcOutIdx)){
      outQcs <- coVar[which(outIdx)[qcOutIdx], 1]
      outQcs <- setdiff(outQcs, emailRec)
      if(length(outQcs) > 0 & !is.null(emailAddress)){
        body <- paste0('The following quality control sample(s) were outlier(s) in the PCA scores:\n',
                       paste0(outQcs, collapse = '\n'), '\nplease check your LC-MS run.')
        # subject
        subject <- "simExTargId warning: check your LC-MS run - QC outlier(s)."
        # add angle brackets to email Address
        emailAddressTmp <- paste0('<', emailAddress, '>')
        email <- lapply(1:length(emailAddressTmp), function(x){
          emailTmp <- try(sendmailR::sendmail(from=emailAddressTmp[1],
                                              to=emailAddressTmp[x],
                                              subject=subject, msg=body,
                                              control=mailControl))})
        # if any error caught then send to console
        if(class(email) == 'try-error'){
          warning(attr(email, 'condition')$message)
        }
        emailRec <- c(emailRec, outQcs)
      }
      }
    # write outlier removed peak table to file
    xcmsPeaks <- pcaOutResults$outRem
    write.csv(outliers, paste0(tmpPcaResultsDir, "/",
                                length(obsNames) - length(remNames),
                                "_outliers_", nSamplesTmp,
                                "_files.csv"), row.names=FALSE)
    } else {
    message("No outliers were removed therefore no outliers table will be saved in the 03.PCA/",
            basename(tmpPcaResultsDir), " sub-directory...\n")
    flush.console()
    }
    # subset covariates table
    remNamesCoVar <- coVar[, 1] %in% remNames
    # subset Covariates table
    coVarOutRem <- coVar[remNamesCoVar, , drop=FALSE]
    # sort peak table in to acquisition order if neccessary
    varCols <- setdiff(colnames(xcmsPeaks), remNames)
    acqOrder <- match(remNames, coVarOutRem[, 1])
    sampClTmp <- rep('unknown', length(acqOrder))
    sampClTmp[!is.na(acqOrder)] <- coVarOutRem[acqOrder[!is.na(acqOrder)], 2]
    # if missing
    acqOrder <- zoo::na.spline(acqOrder)
    remNames <- remNames[order(acqOrder)]
    sampClTmp <- sampClTmp[order(acqOrder)]
    qcNames <- remNames[grepl('^QC$', sampClTmp)]
    xcmsPeaks <- xcmsPeaks[, c(varCols, remNames)]
    # subset the original raw xcms peaks for peak monitoring
    rawColNames <- intersect(c(varCols, remNames), colnames(xcmsPeaksRaw))
    eicIdx <- match(xcmsPeaksRaw$EICno, xcmsPeaks$EICno)
    xcmsPeaksRaw <- xcmsPeaksRaw[!is.na(eicIdx), rawColNames, drop=FALSE]
    # check peak monitoring
    if(!is.null(metab) & length(qcNames) >= 3){
      if(is.data.frame(metab)){
      idxTmp <- match(remNames, coVar[, 1])
      colsPeakMon <- rep('darkgrey', length(idxTmp))
      names(colsPeakMon)[is.na(idxTmp)] <- 'unknown'
      names(colsPeakMon)[!is.na(idxTmp)] <- coVar[idxTmp[!is.na(idxTmp)], 2]
      labPeakMon <- names(colsPeakMon)
      colsPeakMon[grepl('^sample$', names(colsPeakMon))] <- 'black'
      colsPeakMon[grepl('^blank$', names(colsPeakMon))] <- 'blue'
      colsPeakMon[grepl('^QC$', names(colsPeakMon))] <- 'red'
      colsPeakMon[grepl('^ccQC$|^cc$', names(colsPeakMon))] <- 'tomato'
      remGroups <- unique(names(colsPeakMon[colsPeakMon == 'darkgrey']))
      remGroups <- remGroups[remGroups != 'unknown']
      if(length(remGroups) > 0){
        remColIdx <- match(names(colsPeakMon), remGroups)
        colsPeakMon[!is.na(remColIdx)] <- gg_color_hue(4 + length(remGroups))[remColIdx[!is.na(remColIdx)] + 4]
      }
    peakMonitorArgs$launchApp <- FALSE
    forbidIdx <- names(peakMonitorArgs) %in% c("ppm", "rtdev", "maxSignalAtt", "percBelow", 'launchApp')
    peakMonitorArgs <- peakMonitorArgs[which(forbidIdx)]
    peakMonitorArgs$xcmsOutput <- xcmsPeaksRaw
    peakMonitorArgs$obsNames <- remNames
    peakMonitorArgs$qcNames <- qcNames
    # check if any duplicate names in metabolite table
    if(any(dup <- duplicated(metab[, 1]))){
      for(dupname in unique(metab[dup, 1])) {
        dupidx <- which(metab[, 1] == dupname)
        metab[dupidx, 1] <- paste(metab[dupidx, 1], seq(along = dupidx), sep = "__")
      }
    }
    peakMonitorArgs$metab <- metab
    peakMonitorArgs$obsTypeStr <- labPeakMon
    peakMonitorArgs$obsTypeCol <- colsPeakMon

    monitPeaks <- suppressWarnings(do.call(peakMonitor, peakMonitorArgs))
    # if necessary send an email
      if(monitPeaks$warnMessage != ''){
        outQcs <- setdiff(qcNames[length(qcNames)], emailRec)
        if(length(outQcs) > 0 & !is.null(emailAddress)){

          body <- paste0(monitPeaks$warnMessage, '\n\n',
                         paste0(as.vector(t(monitPeaks$summaryTable)), collapse='\n'))
          # subject
          subject <- "simExTargId warning: check your LC-MS run - QC signal attenuated."
          # add angle brackets to email Address
          emailAddressTmp <- paste0('<', emailAddress, '>')
          email <- lapply(1:length(emailAddressTmp), function(x){
            emailTmp <- try(sendmailR::sendmail(from=emailAddressTmp[1],
                                                to=emailAddressTmp[x],
                                                subject=subject, msg=body,
                                                control=mailControl))})
          # if any error caught then send to console
          if(class(email) == 'try-error'){
            warning(attr(email, 'condition')$message)
          }
          emailRec <- c(emailRec, outQcs)
        }
      }
    # save to directory
    peakMonDir <- paste0(outputDir, '/peakMonitor')
    suppressWarnings(dir.create(peakMonDir))
    writeCsvAttempt(df=monitPeaks$monitMetab, fileName=paste0(peakMonDir,
                                                         "/monitorMetaboliteResults.csv"))
    writeCsvAttempt(df=monitPeaks$summaryTable, fileName=paste0(peakMonDir,
                                                         "/summaryTable_", nSamplesTmp, "files.csv"))
     }
    } # end peak-monitor
    # stats analysis
    # 1. reduce number of features by weighted mean of CAMERA pseudospectra
    sumPcGroup <- table(xcmsPeaks$pcgroup)
    message(length(sumPcGroup), ' CAMERA pseudospectra identified.\n')
    flush.console()
    # summed peak area for weighted mean calculation
    xcmsPeaks$peakSum <- apply(LogTransBase ^ xcmsPeaks[, remNames, drop=FALSE], 1, sum)
    # which above min pseudospectrum number
    sumPcGroup <- sumPcGroup[sumPcGroup >= minFeatPseudo]
    message(length(sumPcGroup), ' CAMERA pseudospectra consisted of at least n=',
            minFeatPseudo, ' LC-MS features (minFeatPseudo argument)\n')
    flush.console()

    deconvDf <- xcmsPeaks[xcmsPeaks$pcgroup %in% names(sumPcGroup), , drop=FALSE]
    message('Calculating weighted mean for ', length(sumPcGroup), ' CAMERA pseudospectra.\n')
    flush.console()
    obsIndx <- which(colnames(deconvDf) %in% remNames)

    deconvDf <- by(deconvDf, deconvDf$pcgroup,
                     function(x){
                       wts.tmp <- x[, 'peakSum']
                       wt.mean.tmp <- apply(x[, obsIndx], 2, function(y){
                         weighted.mean(LogTransBase ^ y, wts.tmp)
                       })
                       maxFeat.indx <- which.max(wts.tmp)
                       reqCols <- 1:ncol(x)
                       reqCols <- setdiff(reqCols, obsIndx)
                       nFeatPspec <- nrow(x)
                       names(nFeatPspec) <- "nFeatPspec"
                       wt.mean.tmp <- c(x[maxFeat.indx, reqCols], nFeatPspec,
                                        wt.mean.tmp)
                       return(wt.mean.tmp)
                     })
    deconvDf <- data.frame(do.call(rbind, deconvDf), stringsAsFactors = FALSE)
    deconvDf <- data.frame(lapply(deconvDf, as.character),
                             stringsAsFactors = FALSE)
    # # remove any remaining isotopes
    # message('Removing any remaining isotopes from weighted mean table...\n')
    # flush.console()
    # deconfDf <- deconvDf[grepl('M\\+') == FALSE, , drop=FALSE]
    # make EICno integer
    deconvDf[, 1] <- as.integer(deconvDf[, 1])
    # log transform
    deconvDf[, remNames] <- apply(deconvDf[, remNames], 2, as.numeric)
    deconvDf <- MetMSLine::zeroFill(deconvDf, remNames, zeroFillValue)
    deconvDf <- MetMSLine::logTrans(deconvDf, remNames, LogTransBase)
    # save results
    message("saving weighted mean CAMERA pseudospectra peak table to 02.preProc directory...\n")
    flush.console()
    writeCsvAttempt(df=deconvDf, fileName=paste0(preProcDir,
                                                                "/weightMean_CAMERApseudo_", nSamplesTmp, "files.csv"))

    statRes <- vector('list', ncol(coVarOutRem) - 2)
    statResWMean <- vector('list', ncol(coVarOutRem) - 2)
    statResBatchAdj <- vector('list', ncol(coVarOutRem) - 2)
    statResBatchAdjWMean <- vector('list', ncol(coVarOutRem) - 2)

    message("05. statistical analysis, ", length(statRes), " co-variates...\n")
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
      coVarOutRem <- coVarOutRem[, missingIndx == FALSE, drop=FALSE]
    }
    # create logical data frame of missing values
    coVarOutRem_logi <- data.frame(matrix(TRUE, nrow=nrow(coVarOutRem),
                                          ncol=ncol(coVarOutRem)))
    coVarOutRem_logi[is.na(coVarOutRem)] <- FALSE
    coVarOutRem_logi[coVarOutRem == ""] <- FALSE

    # if nCores then run in parallel
    if(!is.null(nCores)){
      # covariate name
      message("co-variate table columns:\n",
              paste0(colnames(coVarOutRem)[3:ncol(coVarOutRem)], "\n"), "\n",
              nrow(xcmsPeaks), " multiple comparisons pre-processed peak table\n")
      flush.console()

      # start cluster
      message(paste0("Starting SNOW cluster with ", nCores, " local sockets...\n"))
      flush.console()
      cl <- parallel::makeCluster(nCores)
      doSNOW::registerDoSNOW(cl)

      # foreach and dopar from foreach package
      statRes <- foreach(i=3:ncol(coVarOutRem)) %dopar% {
        MetMSLine::coVarTypeStat(peakTable=xcmsPeaks,
                                 obsNames=coVarOutRem[coVarOutRem_logi[, i], 1],
                                 coVariate=coVarOutRem[coVarOutRem_logi[, i], i],
                                 base=LogTransBase,
                                 MTC=pAdjustMethod)
      }
      # stop SNOW cluster
      parallel::stopCluster(cl)

      # covariate name
      message("co-variate table columns:\n",
              paste0(colnames(coVarOutRem)[3:ncol(coVarOutRem)], "\n"), "\n",
              nrow(xcmsPeaks), " multiple comparisons CAMERA pseudospectra weighted mean table\n")
      flush.console()

      # start cluster
      message(paste0("Starting SNOW cluster with ", nCores, " local sockets...\n"))
      flush.console()
      cl <- parallel::makeCluster(nCores)
      doSNOW::registerDoSNOW(cl)

      # foreach and dopar from foreach package
      statResWMean <- foreach(i=3:ncol(coVarOutRem)) %dopar% {
        MetMSLine::coVarTypeStat(peakTable=deconvDf,
                                 obsNames=coVarOutRem[coVarOutRem_logi[, i], 1],
                                 coVariate=coVarOutRem[coVarOutRem_logi[, i], i],
                                 base=LogTransBase,
                                 MTC=pAdjustMethod)
      }
      # stop SNOW cluster
      parallel::stopCluster(cl)
      # if batch adjustment was necessary then stats analysis on other data frame
      if(length(batchAdjusted) > 0){
        message("Performing statistical analyses with batch adjusted data...\n")
        # start cluster
        message(paste0("Starting SNOW cluster with ", nCores, " local sockets...\n"))
        flush.console()
        cl <- parallel::makeCluster(nCores)
        doSNOW::registerDoSNOW(cl)

        # foreach and dopar from foreach package
        statResBatchAdj <- foreach(i=3:ncol(coVarOutRem)) %dopar% {
          MetMSLine::coVarTypeStat(peakTable=batchAdjusted,
                                   obsNames=coVarOutRem[coVarOutRem_logi[, i], 1],
                                   coVariate=coVarOutRem[coVarOutRem_logi[, i], i],
                                   base=LogTransBase,
                                   MTC=pAdjustMethod)
        }
        # stop SNOW cluster
        parallel::stopCluster(cl)

        # covariate name
        message("co-variate table columns:\n",
                paste0(colnames(coVarOutRem)[3:ncol(coVarOutRem)], "\n"), "\n",
                nrow(xcmsPeaks), " multiple comparisons CAMERA pseudospectra weighted mean table\n")
        flush.console()

        # start cluster
        message(paste0("Starting SNOW cluster with ", nCores, " local sockets...\n"))
        flush.console()
        cl <- parallel::makeCluster(nCores)
        doSNOW::registerDoSNOW(cl)

        # foreach and dopar from foreach package
        statResBatchAdjWMean <- foreach(i=3:ncol(coVarOutRem)) %dopar% {
          MetMSLine::coVarTypeStat(peakTable=deconvDf,
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
    pb <- txtProgressBar(min=0, max=length(statRes), style=3)

    for(i in 3:ncol(coVarOutRem)){
      # set progress bar
      setTxtProgressBar(pb, i)
      message("co-variate table column: ", colnames(coVarOutRem)[i], ", ",
              nrow(xcmsPeaks), " multiple comparisons")
      flush.console()
      # auto select coVariate type
      statRes[[i - 2]] <- MetMSLine::coVarTypeStat(xcmsPeaks,
                                                   obsNames=coVarOutRem[coVarOutRem_logi[, i], 1],
                                                   coVariate=coVarOutRem[coVarOutRem_logi[, i], i],
                                                   base=LogTransBase,
                                                   MTC=pAdjustMethod)
      statResWMean[[i - 2]] <- MetMSLine::coVarTypeStat(deconvDf,
                                                        obsNames=coVarOutRem[coVarOutRem_logi[, i], 1],
                                                        coVariate=coVarOutRem[coVarOutRem_logi[, i], i],
                                                        base=LogTransBase,
                                                        MTC=pAdjustMethod)
      # if necessary then stats analysis batchadjusted data
      if(length(batchAdjusted) > 0){
        message("Batch Adjusted, co-variate table column: ", colnames(coVarOutRem)[i], ", ",
                nrow(xcmsPeaks), " multiple comparisons")
        flush.console()

      statResBatchAdj[[i - 2]] <- MetMSLine::coVarTypeStat(batchAdjusted,
                                                           obsNames=coVarOutRem[coVarOutRem_logi[, i], 1],
                                                           coVariate=coVarOutRem[coVarOutRem_logi[, i], i],
                                                           base=LogTransBase,
                                                           MTC=pAdjustMethod)
      statResBatchAdjWMean[[i - 2]] <- MetMSLine::coVarTypeStat(batchAdjusted,
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
    save(statRes, file=paste0(rDataDir, "/statResults.RData"))
    save(statResWMean, file=paste0(rDataDir, "/statResults_weightMean.RData"))

    # if necessary save batch adjusted statistical test results
    if(length(batchAdjusted) > 0){
      message('saving batch adjusted statistical analysis objects (.RData) in RDataFiles directory...')
      flush.console()
      save(statResBatchAdj, file= paste0(rDataDir, "/statResultsBatchAdj.RData"))
      save(statResBatchAdjWMean, file= paste0(rDataDir, "/statResultsBatchAdj_weightMean.RData"))
    }

    # covar names
    coVarNames <- colnames(coVarOutRem)[3:ncol(coVarOutRem)]
    # add stat results output to peak table
    statRes <- do.call(cbind, lapply(1:length(statRes), function(x){
          statResTmp <- statRes[[x]]
          if(!is.character(statResTmp)){
          colNamesTmp <- paste(statResTmp$method, statResTmp$MTC, coVarNames[x],
                               sep="_")
          resDfTmp <- statResTmp$result
          colnames(resDfTmp) <- paste0(colnames(resDfTmp), "_", colNamesTmp)
          return(resDfTmp)
          }}))
    statResWMean <- do.call(cbind, lapply(1:length(statResWMean), function(x){
      statResTmp <- statResWMean[[x]]
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
    statRes <- data.frame(xcmsPeaks[, remColsIndx], statRes,
                              xcmsPeaks[, obsIndx], stringsAsFactors=FALSE)
    # if necessary extract batch adjusted test results
    if(length(batchAdjusted) > 0){
      statResBatchAdj <- do.call(cbind, lapply(1:length(statResBatchAdj),
                                                   function(x){
          statResTmp <- statResBatchAdj[[x]]
          if(!is.character(statResTmp)){
          colNamesTmp <- paste("batchAdj", statResTmp$method, statResTmp$MTC, coVarNames[x],
                               sep="_")
          resDfTmp <- statResTmp$result
          colnames(resDfTmp) <- paste0(colnames(resDfTmp), "_", colNamesTmp)
          return(resDfTmp)
          }}))
      # insert stats res and write to results directory
      statResBatchAdj <- data.frame(xcmsPeaks[, remColsIndx],
                                        statResBatchAdj, xcmsPeaks[, obsIndx],
                                        stringsAsFactors=FALSE)
    } # cond batch adjusted

    # insert stats analysis into peak table
    obsIndx <- match(coVarOutRem[, 1], colnames(deconvDf))
    # remaining columns
    remColsIndx <- setdiff(1:ncol(deconvDf), obsIndx)
    # insert stats res and write to results directory
    statResWMean <- data.frame(deconvDf[, remColsIndx], statResWMean,
                               deconvDf[, obsIndx], stringsAsFactors=FALSE)

    if(length(batchAdjusted) > 0){
      statResBatchAdjWMean <- do.call(cbind, lapply(1:length(statResBatchAdjWMean),
                                               function(x){
                                                 statResTmp <- statResBatchAdjWMean[[x]]
                                                 if(!is.character(statResTmp)){
                                                   colNamesTmp <- paste("batchAdj", statResTmp$method, statResTmp$MTC, coVarNames[x],
                                                                        sep="_")
                                                   resDfTmp <- statResTmp$result
                                                   colnames(resDfTmp) <- paste0(colnames(resDfTmp), "_", colNamesTmp)
                                                   return(resDfTmp)
                                                 }}))
      # insert stats res and write to results directory
      statResBatchAdjWMean <- data.frame(deconvDf[, remColsIndx],
                                         statResBatchAdj, deconvDf[, obsIndx],
                                         stringsAsFactors=FALSE)
      } # cond batch adjusted


     message("Writing statistical results output, ", nSamplesTmp, " samples...\n")
     flush.console()
     # write csv results
     writeCsvAttempt(df=statRes, fileName=paste0(statsDir, "/statResults.csv"))
     writeCsvAttempt(df=statResWMean, fileName=paste0(statsDir, "/statResults_wMean.csv"))
     message("...done\n")
     flush.console()

    if(length(batchAdjusted) > 0){

    message("Writing batch adjusted statistical results output, ", nSamplesTmp, " files...\n")
    flush.console()
    # write csv results
    writeCsvAttempt(df=statResBatchAdj, fileName=paste0(statsDir,
                                                         "/statResults_BatchAdj.csv"))
    writeCsvAttempt(df=statResBatchAdjWMean, fileName=paste0(statsDir, "/statResults_BatchAdj_wMean.csv"))

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
  # if(all(SampleConv == T)){
  #
  # message("Generating final diffreport and EICs...\n")
  # flush.console()
  # setwd(rDataDir)
  # # generate diff report
  # suppressWarnings(do.call(xcms::diffreport, xcmsDiffRepArgs))
  #
  # message("finished...\n")
  # flush.console()
  # break
  # }
} # end while loop
} else {
  message("all of the samples have already been converted to mzXML files stopping")
  flush.console()
  }# end if all samples already converted
} # end function
