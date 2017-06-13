## ---- include=F----------------------------------------------------------
library(simExTargId)

## ---- eval=FALSE---------------------------------------------------------
#  # select directory where example directories for raw files can be created and also results
#  # output saved
#  # For example the C: directory
#  # create a dummy raw MS1 profiling data directory
#  # for example
#  studyName <- paste0(gsub(" ", "", format(Sys.time(), "%Y %m %d")), "_simExTargId_Example")
#  dummyRawDir <- paste0("C:/", studyName)
#  # create directory
#  dir.create(dummyRawDir)
#  # replace email address with a vector of one or more email addresses of
#  # people in your laboratory group
#  emailNotifier(rawDir=dummyRawDir, emailAddress='johndoe@emailprovider.com',
#                emailTime=5)

## ---- eval=FALSE---------------------------------------------------------
#  # simExTargId extdata directory
#  extdataDir <- system.file("extdata", package="simExTargId")
#  # list blank files
#  blanksRaw <- list.files(extdataDir, pattern="blank", full.names=TRUE)
#  # list plasma IPA extract files
#  samplesRaw <- list.files(extdataDir, pattern="sample", full.names=TRUE)
#  # covariates file
#  coVariates <- paste(extdataDir, "coVariates.csv", sep="/")
#  # illustrative metabolite database table
#  metabDb <- paste(extdataDir, 'exampleMetabDatabase.csv', sep='/')
#  # identify number of virtual cores for parallel processing using parallel package
#  nCores <- parallel::detectCores()
#  
#  # as no qc files then use sample files twice to illustrate the
#  # peakMonitor function
#  
#  # move the blank files into the temporary directory to start the process
#  blankRawCopies <- paste(dummyRawDir, basename(blanksRaw), sep="/")
#  file.copy(from=blanksRaw, to=blankRawCopies)
#  # set the file time to simulate the files having been acquired at least 5 mins
#  # since last modification
#  setTheTime <- function(fileCopy, time){
#    Sys.setFileTime(fileCopy, Sys.time() - time)
#  }
#  # apply to newly copied files
#  sapply(blankRawCopies, setTheTime, 240)
#  
#  # move the plasma samples twice first time rename as QC
#  # and set the file time less than 5 mins
#  samplesRawCopies <- paste(dummyRawDir, basename(samplesRaw), sep="/")
#  file.copy(from=samplesRaw, to=samplesRawCopies)
#  qcFiles <- gsub('sample', 'qc', samplesRawCopies)
#  file.rename(from=samplesRawCopies, to=qcFiles)
#  file.copy(from=samplesRaw, to=samplesRawCopies)
#  
#  # apply to newly copied files
#  sapply(c(samplesRawCopies, qcFiles), setTheTime, 300)
#  
#  # Start simExTargId function
#  simExTargId(rawDir=dummyRawDir, studyName = studyName, analysisDir='C:/',
#              coVar=coVariates, metab=metabDb, nCores=nCores, ionMode='nega',
#              minFiles=3)

## ---- eval=FALSE---------------------------------------------------------
#  peakMonitor(analysisDir=paste0(dummyRawDir, "_analysis/NEG/output/peakMonitor"))

## ----eval=FALSE----------------------------------------------------------
#  # this command will open the application in your web-browser
#  targetId(analysisDir=paste0(dummyRawDir, "_analysis/NEG/output/04.stats"))

