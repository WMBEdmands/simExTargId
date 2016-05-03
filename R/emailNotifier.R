#' email notification system for LC-MS analyses
#' @param rawDir character full path name of raw data directory into which raw data will 
#' eventually/ has already been started to be written. 
#' @param emailTime numeric length of time (in minutes) after the last raw data 
#' file was written before an email notification will be sent. minimum value is 5 minutes.
#' @param mailControl List of SMTP server settings see \code{\link{sendmail}} for details. 
#' Example given is for google mail.
#' @param emailAddress character email address from and to which to send warning email
#' that run may have stopped. (if not supplied then email notification will not be sent) 
#' see \code{\link{sendmail}}.
#' 
#' @export
emailNotifier <- function(rawDir=NULL, emailAddress=NULL, emailTime=20, 
                          maxTime=60, 
                          mailControl=list(smtpServer="ASPMX.L.GOOGLE.COM")){
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
  rawFiles <- dir(rawDir, full.names=T, pattern="\\.d$|\\.RAW$")
  
  
  if(length(rawFiles) == 0){
    # wait for first file to be acquired
    message("waiting for the first raw data file to be acquired...\n")
    flush.console()
  }
  while(length(rawFiles) == 0){
    # detect raw samples in directory
    rawFiles <- dir(rawDir, full.names=T, pattern="\\.d$|\\.RAW$")
    if(length(rawFiles) == 1){
      break
    }
  }
  # email sent counter to prevent more than one email being sent.
  emailSent <- 0
  
  while(length(rawFiles) > 0){
    # put system to sleep for 5 minutes before checking raw file status
    Sys.sleep(300)
    # detect raw samples in directory
    rawFiles <- dir(rawDir, full.names=T, pattern="\\.d$|\\.RAW$")
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
    
    # test to see when last raw file was modified and then send email if necessary
    if(min(timeChIndx) > emailTime){
      if(emailSent == 0){
        # iterate email sent
        emailSent <- 1
        # if email address supplied then send a warning to that email, from that email
          # email body
          body <- paste0('The last raw data file ', basename(names(which.min(timeChIndx))),
                         ' was modified ',  round(min(timeChIndx), digits=2), 
                         ' minutes ago, please check your LC-MS run.')
          # subject
          subject <- "simExTargId warning: check your LC-MS run!"
          # add angle brackets to email Address
          emailAddress <- paste0('<', emailAddress, '>')
          email <- lapply(1:length(emailAddress), function(x){
            emailTmp <- try(sendmailR::sendmail(from=emailAddress[1], 
                                                to=emailAddress[x], 
                                                subject=subject, msg=body, 
                                                control=mailControl))})
          # if any error caught then send to console
          if(class(email) == 'try-error'){
            warning(attr(email, 'condition')$message)  
          }
       
      }
  }
    # if longer than an hour then stop the process.
    if(min(timeChIndx) > maxTime){
      stop("The last raw data file was written longer than an hour ago...stopping the email notification process")
    }
  }
}