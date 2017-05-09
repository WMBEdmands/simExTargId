#' Convert proprietory instrument manufacturer file types to .mzXML files
#' 
#' @details MSConvert (ProteoWizard \url{http://proteowizard.sourceforge.net/downloads.shtml}) 
#' must be installed and in the case of windows
#' in the path in the environmental variables. The mzXML file conversion will occur
#' through shell commands using the command line version of MSConvert. If more than
#' one cpu is utilized the mzXML file conversion process will occur as a parallel
#' computation.
#' @param rawFiles character vector of full path names to raw data files.
#' @param mzXmlDir character full path to directory into which to write the mzXml
#' files output by MSConvert.
#' @param nCores numeric number of computer cores for parallel computation.
#' @param subSetSecs numeric vector of a minimum and maximum time window in 
#' mzXML files (e.g. c(1500, 1900)).
#' @param centroid do the raw data files need to be centroided during conversion by MSConvert
#' \url{http://proteowizard.sourceforge.net/downloads.shtml}. NB. centroiding of data
#' is necessary for utilization of the "centWave" algorithm of xcms (\code{\link{findPeaks.centWave-methods}}).
#' @param zlib logical whether or not to apply zlib compression (default = TRUE). Will reduce
#' file size but may affect downstream data analysis.
#' @param MS2 logical should raw LC-MS data files converted by MSConvert to mzXml be converted
#' as MS/MS files (e.g. data-dependent MS/MS files).
#' @export
rawFileConvert <- function(rawFiles=NULL, mzXmlDir=NULL, nCores=NULL, subSetSecs=NULL,
                           zlib=T, centroid=TRUE, MS2=FALSE){
  #error handling
 if(is.null(rawFiles)){
   stop('Argument rawFiles is missing with no default')
 }
  # if necessary select analysis (i.e. results output) directory 
  if(is.null(mzXmlDir)){
    message("Select the directory your mzXML files will be written to...\n")
    flush.console()
    
    mzXmlDir <- tcltk::tk_choose.dir(default = "", 
                                        caption = "Select the directory your .mzXML files/ results files will be written to...")
  }
  if(!is.null(nCores)){
    if(!require(foreach)){
      stop('foreach package must be installed. see ?install.packages')
    }
  }
  # Establish msconvert commands to be sent to shell commands, centroid/ MS2 
  # (i.e. data dependent/independent) if necessary
  if(centroid == T){
    convType <- ifelse(MS2 == T, paste0('" --32 --mzXML', ifelse(zlib == T, ' --zlib ', ' '), 
                       '--filter "peakPicking true 1-2"'),
                       paste0('" --32 --mzXML', ifelse(zlib == T, ' --zlib ', ' '), 
                       '--filter "peakPicking true 1-"'))
  } else {
    convType <- ifelse(MS2 == T, paste0('" --32 --mzXML', ifelse(zlib == T, ' --zlib ', ' '), 
                       '--filter "msLevel 1-2"'),
                       paste0('" --32 --mzXML', ifelse(zlib == T, ' --zlib ', ' '), 
                       '--filter "msLevel 1-"'))
  }
  # if necessary subset by min and max time in seconds
  if(!is.null(subSetSecs)){
    # error handling, if length of argument is not equal two (i.e. min and max)
    if(length(subSetSecs) != 2){
      stop("subSetSecs argument must be a numeric vector of length two")
    } else {
    convType <- paste0(convType, ' --filter "scanTime [', subSetSecs[1], ',', subSetSecs[2], ']"') 
    }
  }
  
  message("Converting to mzXML and initial peak picking: \n", 
          paste(basename(rawFiles), collapse="\n"), "\n") 
  flush.console()
  
  command <- paste0('msconvert "', rawFiles, convType, ' -o "', mzXmlDir, '"')
  # if nCores not null then convert to mzXML in parallel
  if(!is.null(nCores) & length(command) > 1){
    # start cluster
    message(paste0("Starting SNOW cluster with ", nCores, " local sockets...\n"))
    flush.console()
    cl <- parallel::makeCluster(nCores) 
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
}