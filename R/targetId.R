#' Simultaneous MS1 profiling and target identification visualization
#'
#' @description Visualization of deconvoluted data obtained during an MS1-level
#' profiling experiment. All statResults (.csv) files within an analysis directory
#' will be read and can be visualized as a volcano plot and as reactive
#' tables indicating the direction of the median fold change. This allows the
#' experimenter to visualize the emergence of statistically significant
#' pseudospectra features before the completion of a MS1-level profiling experiment.
#' minimum fold change and maximum multiple testing adjusted p-/q-values (\code{\link{coVarStatType}}) can be
#' set by the experimenter. When a list of possible targets has been determined
#' the table can be downloaded/ saved as a Csv file. This list of statistically
#' significant targets can therefore be readily included in a MS/MS fragmentation
#' experiment before the end of a profiling experiments completion.
#'
#' @param analysisDir (i.e. 04.stats)the full path of the analysis directory containing the
#'  (.csv) files output from \code{\link{simExTargId}}.
#'
#' @seealso \code{\link{simExTargId}}.
#'
#' @export
targetId <- function(analysisDir=NULL){

  if(is.null(analysisDir)){
    message("window opened to select analysis (/04.stats) results directory containing the results (.csv) files...")
    flush.console()

    analysisDir <- tcltk::tk_choose.dir(default = "",
                                        caption = "Select the analysis (/04.stats) results directory containing the results (.csv) files...")
  }
  if(basename(analysisDir) != '04.stats'){
    stop('The analysisDir must be names "04.stats".\n')
  }
  analysisDir <<- analysisDir
  appDir <- system.file("shiny-apps", "targetId", package = "simExTargId")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `simExTargId`.", call. = FALSE)
  }
  source(paste0(appDir, '/global.R'))
  source(paste0(appDir, '/ui.R'))
  source(paste0(appDir, '/server.R'))

  shiny::runApp(list(ui=targIdui, server=targIdserver), launch.browser=TRUE)

} # end function
