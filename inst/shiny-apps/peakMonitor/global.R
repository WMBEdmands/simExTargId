# calculate PCA plot of raw data and plot using rCharts
if(!require(pcaMethods)){
  stop('The package pcaMethods must be installed to use the peakMonitor application...')
}
if(!require(data.table)){
  stop('The package data.table must be installed to use the peakMonitor application...')
}
if(!require(DT)){
  stop('The package DT must be installed to use the peakMonitor application...')
}
if(!require(MetMSLine)){
  stop('The package MetMSLine must be installed to use the peakMonitor application...')
}

exMetabFilePath <- paste0(analysisDir, '/monitorMetaboliteResults.csv')
if(!file.exists(exMetabFilePath)){
  stop('The expected file:\n', exMetabFilePath, '\nWas not found in the analysis directory. This file should not be renamed and no columns modified to be read by peakMonitor.\n')
}
exMetab <- as.data.frame(data.table::fread(exMetabFilePath, stringsAsFactors = FALSE, header = TRUE))
# observation names
obsNames <- colnames(exMetab)[{grep('belowMedianSigAtt', colnames(exMetab)) + 1}:ncol(exMetab)]
obsTypeStr <- unlist(exMetab[1, obsNames])
obsTypeCol <- unlist(exMetab[2, obsNames])
exMetab <- exMetab[-2:-1, , drop=FALSE]
exMetab[, 1] <- as.integer(exMetab[, 1])
exMetab[, obsNames] <- apply(exMetab[, obsNames], 2, as.numeric)
exMetabLog <- MetMSLine::zeroFill(exMetab, obsNames)
exMetabLog <- MetMSLine::logTrans(exMetabLog, obsNames)
pcaRes <- pcaMethods::pca(t(exMetabLog[, obsNames]), nPcs = 2, method="svd",
                          scale="pareto", centre=TRUE, cv="q2")
# Get the estimated scores
scoresDf <- as.data.frame(pcaRes@scores)
pcaSumm <- invisible(summary(pcaRes))
# Get the loadings
loadingsDf <- as.data.frame(pcaRes@loadings)

# add number padded ids
paddSampTypes <- obsTypeStr
for(nameIt in unique(obsTypeStr)){
  dupidx <- which(obsTypeStr == nameIt)
  paddSampTypes[dupidx] <- paste(obsTypeStr[dupidx], sprintf('%04d', seq(along = dupidx)), sep = "_")
}

qcNames <- obsNames[grep('^QC$', obsTypeStr)]
# rearrange column order
ordCols <- c(1:4, grep('obs_mzmed', colnames(exMetab)):ncol(exMetab))
ordCols <- c(ordCols, setdiff(1:ncol(exMetab), ordCols))
exMetab <- exMetab[, ordCols]
