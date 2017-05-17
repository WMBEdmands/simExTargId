# identify results files
resFileNames <- list.files(analysisDir, full.names=TRUE, pattern="\\.csv",
                           recursive=TRUE)
if(length(resFileNames) == 0){
  stop('no statistical analysis result files (.csv) were identified in the directory:\n',
       analysisDir, '\nPlease check and try again.\n')
}
# add all results files into list
message("Reading ", length(resFileNames), " (.csv) results files please wait..")
flush.console()
resFiles <- lapply(resFileNames, function(x) as.data.frame(data.table::fread(x, header=TRUE, stringsAsFactors=FALSE)))
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
