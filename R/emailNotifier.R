#' email notification system for LC-MS analyses
#' @param rawDir character full path name of raw data directory into which raw data will
#' eventually/ has already been started to be written.
#' @param emailAddress character vector email address(es) from and to which to send warning email
#' that run may have stopped. (if not supplied then email notification will not be sent)
#' see \code{\link{sendmail}}. The first email address in the character vector will
#' be the email address which warning emails are sent from to all of the other
#' email addresses supplied. N.B. Currently the email addresses must use the same SMTP
#' server setting (e.g. all google mail, see mailControl argument below).
#' @param emailTime numeric length of time (in minutes) after the last raw data
#' file was written before an email notification will be sent. minimum value is 5 minutes.
#' @param maxTime numeric the maximum time before the email notifier process will stop
#' altogether so it does not continue ad infinitum and the eventual return of
#' the cosmos to a state of maximum entropy (default = 60 mins).
#' @param mailControl List of SMTP server settings see \code{\link{sendmail}} for details.
#' Example given is for google mail.
#' @param minFiles integer the minimum number of files that must be collected
#' before the median absolute deviation is calculated (default = 5).
#' @param nMad numeric the number of median absolute deviations smaller a file
#' must be for an email warning to be sent (default = 3). Reducing this
#' value may make it more sensitive to subtle reductions in file size but increase
#' the risk of false positive warning email messages. Utilizes the base function
#' \link{mad} in the stats package.
#'
#' @export
emailNotifier <- function(rawDir=NULL, emailAddress=NULL, emailTime=10,
                          maxTime=60,
                          mailControl=list(smtpServer="ASPMX.L.GOOGLE.COM"),
                          minFiles=5, nMad=3){
  # error handling
  if(is.null(emailAddress)){
    stop('emailAddress argument is missing with no default')
  } else if(emailTime < 5){
    stop('minimum allowable emailTime value is 5 minutes')
  }

  # if necessary select raw data directory
  if(is.null(rawDir)){
    message("Select the directory your raw data files will be/ have already started to be written to...\n")
    flush.console()

    rawDir <- tcltk::tk_choose.dir(default = "",
                                   caption = "Select the directory your raw data files will be written to...")
  }

  message('monitoring raw file directory for an unexpected stop in the LC-MS acquisition run...')
  flush.console()
  # detect raw samples in directory
  rawFiles <- dir(rawDir, full.names=TRUE, pattern="\\.d$|\\.RAW$")

  # empty vector for datafiles collected more than email time ago in directory
  finishedFiles <- as.character()
  fileSizes <- as.numeric()
  emailSentFiles <- as.character()
  emailSentSize <- as.character()
  if(length(rawFiles) == 0){
    # wait for first file to be acquired
    message("waiting for the first raw data file to be acquired...\n")
    flush.console()
  }
  while(length(rawFiles) == 0){
    # detect raw samples in directory
    rawFiles <- dir(rawDir, full.names=TRUE, pattern="\\.d$|\\.RAW$")
    if(length(rawFiles) == 1){
      break
    }
  }
  # last file collected timer
  lastFile <- 0
  nFilesCompleted <- 0
  while(lastFile < maxTime){
    # put system to sleep for 60 secs before checking raw file status
    Sys.sleep(30)
    # detect raw samples in directory
    rawFiles <- dir(rawDir, full.names=TRUE, pattern="\\.d$|\\.RAW$")
    rawFiles <- setdiff(rawFiles, finishedFiles)
    if(length(rawFiles) > 0){
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
    # file sizes
    fileSizesRaw <- sapply(rawFiles, function(rawFile){
      # if the raw file is a directory then identify all files within
      dirFilesTmp <- list.files(path=rawFile, recursive=TRUE, full.names=TRUE)
      if(length(dirFilesTmp) > 0){
        sizeTmp <- max(file.info(dirFilesTmp)$size)
      } else {
        sizeTmp <- max(file.info(rawFile)$size)
      }
      return(sizeTmp)})
    # if time more than email time then add to finished files vector
    if({nNewFiles <- length(finishedFiles) - nFilesCompleted} > 0){
    cat(paste0(length(finishedFiles), ' (',
               gsub('.+\\.', '.', basename(finishedFiles[1])),
               ') file(s) collected.\n'))
    nFilesCompleted <- nFilesCompleted + nNewFiles

    }
    lastFile <- min(timeChIndx)
    finIdx <- timeChIndx > {emailTime * 1.1}
    finishedFiles <- c(finishedFiles, names(timeChIndx[finIdx]))
    fileSizes <- c(fileSizes, fileSizesRaw[finIdx])
    # test to see when last raw file was modified and then send email if necessary
    if(lastFile > emailTime){
      lastRawFile <- basename(names(which.min(timeChIndx)))
      lastRawFile <- setdiff(lastRawFile, emailSentFiles)
      if(length(lastRawFile) > 0){
        # if email address supplied then send a warning to that email, from that email
          # email body
          body <- paste0('The last raw data file ', lastRawFile,
                         ' was modified ',  round(min(timeChIndx), digits=2),
                         ' minutes ago, please check your LC-MS run.')
          # subject
          subject <- "simExTargId warning: check your LC-MS run!"
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
        emailSentFiles <- c(emailSentFiles, lastRawFile)
      }
    }

    # if file sizes greater than 3 median absolute deviations from median after min 5 files then email
    if(length(fileSizes) > minFiles){
    madFS <- mad(fileSizes)
    medFS <- median(fileSizes)
    loOut <- medFS - {nMad * madFS}
    if(any(outFileSize <- fileSizes <= loOut)){
      outSizeTmp <- basename(finishedFiles[outFileSize])
      outSizeTmp <- setdiff(outSizeTmp, emailSentSize)
      if(length(outSizeTmp) > 0){
        # if email address supplied then send a warning to that email, from that email
        # email body
        body <- paste0('The file size of the following completed raw data file\n ', outSizeTmp,
                       ' is less than ', nMad, ' median absolute deviations smaller than the median file size. please check your LC-MS run.')
        # subject
        subject <- "simExTargId warning: check your LC-MS run (small file size)!"
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
        emailSentSize <- c(emailSentSize, outSizeTmp)
      }
     }
    }
    }
    # if longer than an hour then stop the process.
    if(lastFile > maxTime){
      stop("The last raw data file was written longer than an hour ago...stopping the email notification process")
      break
    }
  }
}# end function
