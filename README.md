simExTargId 
===========
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.806838.svg)](https://doi.org/10.5281/zenodo.806838) Development version 0.2.1 archived on the Zenodo repository. 

Performs simultaneous data conversion, preprocessing and statistical analysis and instrument stoppage/signal drift monitoring throughout the course of a metabolomic profiling experiment. The **simExTargId** package performs raw data to mzXML conversion ([MSConvert](http://proteowizard.sourceforge.net/)), peak-picking, retention time alignment, grouping ([xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html)), ESI adducts/in-source fragment artefact identification ([CAMERA](https://bioconductor.org/packages/release/bioc/html/CAMERA.html)), peak table output pre-processing, automatic outlier detection by PCA and automatic  statistical analysis dependent on co-variates supplied ([MetMSLine](https://github.com/WMBEdmands/MetMSLine)). Signal drift monitoring and possible MS2 target list identification is facilitated during a metabolomic profiling experiment by two shiny user interfaces. An email notification system using the [sendmailR](https://cran.r-project.org/web/packages/sendmailR/index.html)
package is implemented in *simExTargId* to warn the user(s) of instrumental stoppages, outlying quality control samples and signal attenuation/drift. The *simExTargId* package is still in active development but nevertheless should still be at a functional and useful stage for many metabolomic investigators.

Installation
===============
The latest development version and all package dependencies can be installed with one-line of code directly from GitHub using the devtools package. First ensure devtools is installed, instructions can be found here: https://github.com/hadley/devtools
```{r}
devtools::install_github('WMBEdmands/simExTargId', dependencies=c("Depends", "Imports", "Suggests"))
```

Vignette
========
The vignette can be viewed [here](http://bit.ly/2rUQSAk) with example raw data acquired on a Thermo FT-ICR mass spectrometer (this requires playing with the creation time of the raw data files to simulate a real-time data collection).

The *simExTargId* package has thus far only been tested with Agilent (.d) and Thermo (.RAW/.raw) data files on a computer running Windows but depending on interest could be readily extended to other instrument manufacturers and may well already function on an OSX or Linux based operating system.

Details 
=======
The R package utilizes MSConvert ([ProteoWizard](http://proteowizard.sourceforge.net/)), [xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html), [CAMERA](https://bioconductor.org/packages/release/bioc/html/CAMERA.html), [MetMSLine](https://github.com/WMBEdmands/MetMSLine), [Shiny](https://shiny.rstudio.com/) and many other packages to attempt to implement real-time **sim**ultaneous metabolomic MS1-profiling **Ex**periment and statistically relevant MS2 **Targ**et **Id**entification.

A major impetus for development of this package was to provide an email-based early warning system for LC-MS instrumental stoppages/errors but also for more subtle changes such as instrument drift and (PCA-based) outlying pooled quality control samples for example. When collecting a dataset of precious/limited samples such as those with a low volume/quantity provided (e.g. mouse sera, dried blood spot punches) it can be particularly poignant if un-noticed instrumental stoppages leads to degradation of your samples or unwanted variation/batch effects.
Instrumental drifts and potentially MS2 targets can be identified in real time using the two shiny apps developed for *simExTargId* `peakMonitor` and `targetId`.

Narrowing the temporal gap to truly "online" target feature MS2 fragmentation is also a primary goal of *simExTargId*. 

All peak-picking, retention time alignment and grouping is performed by [xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html), then ESI adducts and isotopes are detected by [CAMERA](https://bioconductor.org/packages/release/bioc/html/CAMERA.html) in real time. Following this pre-processing of the peak-table is performed by [MetMSLine](https://github.com/WMBEdmands/MetMSLine). Automatic PCA-based outlier detection is performed and a warning email sent to user(s) if a QC sample is outlying. Furthermore, automatic real-time univariate co-variate based statistical analysis is then performed. 

![example_PCA_Scores_Iteration_1+2](https://github.com/WMBEdmands/simExTargId/blob/master/inst/extdata/pcaOutId_scores_iter1_2.png)
**Example PCA scores iterations 1 and 2** - averaged pooled Quality Control samples should cluster tightly near the origins of the first and second principal components

******
The *simExTargId* workflow can be left running during the MS1-profiling data collection providing a degree of reassurance that serious instrumental difficulties will not go unnoticed (even whilst present in the laboratory and otherwise distracted) and also outlying samples and statistically relevant LC-MS feature targets can be identified and additional experiments/reinjections appended to the experimental worklist before the end.
******

![example_peakMonitorEmail](https://github.com/WMBEdmands/simExTargId/blob/master/inst/extdata/examplePeakMonitorEmail.PNG)

**Example signal attenuation warning email** - pooled quality controls can be used to monitor a database (.csv) of previously identified metabolites if a maximum percentage (>20%) of the monitored metabolites are greater than 20% attenuated (default can be altered by user) then a warning email is sent. The output of the `peakMonitor` function can be visualized in a shiny application.

******
*SimExTargId* also attempts to inculcate a rigorous (and broadly accepted) experimental design and also a directory sub-structure for each experiment to help organize data and files for the user. Sub-directories, tables and plots are then generated in real time. For example, each time a new iteration of *simExTargId* takes place updated PCA plots and tables will appear in each sub-directory. The idea is to preserve the reproducibility by the recording of parameters used (and session information) and providing an easy to navigate and intuitive sub-directory structure.
******

<img align="left" src="https://github.com/WMBEdmands/simExTargId/blob/master/inst/extdata/subdirStr_example.PNG">

**Example sub-directory structure generated by simExTargId upon initiation** - the "data" directory should contain files which should not be modified such as the worklist/co-variates table received from the collaborators or the worklist, the "doc" directory should contain all documentation i.e. manuscript files and method details, the "output" directory contains all output generated by *simExTargId* or by R code created after data collection, the "R" directory is used for all R code and workspace (.RData/.RDS) files. All mzXML files are saved in a separate directory due to their large size and to keep the analysis directory small.
******
A technical note describing the package is currently under-review and awaiting decision until then if you use *simExTargId* please cite the DOI provided by the Zenodo archive.

Licence
=============
The *simExTargId* package is licenced under the GPLv3 (http://www.gnu.org/licenses/gpl.html).

 
