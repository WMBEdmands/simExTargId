#' try to write a csv and pause and alert user if open
#' @param df dataframe
#' @param fileName complete path to new .csv file
writeCsvAttempt <- function(df=NULL, fileName=NULL){
  # write batchAdjusted residual xcms peak table to a .csv
  writeAttempt <- suppressWarnings(try(write.csv(df, fileName, 
                                                 row.names=F), silent=T))
  errorMessage <- 0
  while(class(writeAttempt) == "try-error"){
    if(errorMessage == 0){
      message('Please close the ', basename(fileName), 
              ' file so simExTargId can continue...')
      flush.console()
      errorMessage <- 1
    }
    writeAttempt <- suppressWarnings(try(write.csv(df, fileName, 
                                                   row.names=F), silent=T))
  } 
}