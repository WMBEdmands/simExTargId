# simExTargId (v0.2.0)
Performs simultaneous raw data to mzXML conversion ([MSConvert](http://proteowizard.sourceforge.net/)), peak-picking, retention time alignment, grouping ([xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html)), ESI adducts/in-source fragment artefact identification ([CAMERA](https://bioconductor.org/packages/release/bioc/html/CAMERA.html)), peak table output pre-processing, automatic outlier detection by PCA and automatic  statistical analysis dependent on co-variates supplied ([MetMSLine](https://github.com/WMBEdmands/MetMSLine)), signal drift monitoring and possible MS2 target list identification and real-time shiny-based visualizations are also implemented in real-time during a metabolomic profiling experiment.

The R package utilizes MSConvert ([ProteoWizard](http://proteowizard.sourceforge.net/)), [xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html), [CAMERA](https://bioconductor.org/packages/release/bioc/html/CAMERA.html), [MetMSLine](https://github.com/WMBEdmands/MetMSLine), [Shiny](https://shiny.rstudio.com/) and many other packages to attempt to implement real-time **sim**ultaneous metabolomic MS1-profiling **Ex**perimentation and statistically relevant MS2 **Targ**et **Id**entification.

A main impetus for development of this package was to provide an email-based early warning system for LC-MS instrumental stoppages/errors but also for more subtle changes such as instrument drift and (PCA-based) outlying pooled quality control samples for example. When collecting a dataset of precious/limited samples such as those with a low volume/quantity provided (e.g. mouse sera, dried blood spot punches) it can be particularly poignant if un-noticed instrumental stoppages leads to degradation of your samples or unwanted variation/batch effects.
Instrumental drifts and potentially MS2 targets can be identified in real time using the two shiny apps developed for **simExTargId** *peakMonitor* and *targetId*.

Narrowing the gap to truly "online" target feature MS2 fragmentation is also a primary goal of **simExTargId**. 

All peak-picking, retention time alignment and grouping is performed by [xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html), then ESI adducts and isotopes detected by [CAMERA](https://bioconductor.org/packages/release/bioc/html/CAMERA.html) in real time. Following this pre-processing of the peak-table is performed by [MetMSLine](https://github.com/WMBEdmands/MetMSLine). Automatic PCA-based outlier detection is performed and a warning email sent to user(s) if a QC sample is outlying. Furthermore, automatic real-time univariate co-variate based statistical analysis is then performed. 

The **simExTargId** workflow can be left running during the MS1-profiling data collection providing a degree of reassurance that serious instrumental difficulties will not go unnoticed (whilst not present in the laboratory) and also outlying samples and statistically relevant LC-MS feature target can be identified and additional experiments/reinjected appended to the experimental worklist before the end.

**SimExTargId** also attempts to inculcate a rigorous (and broadly accepted) experimental design and also a organizational directory sub-structure for an experiment, helping to organize data and files for the user. Sub-directories, tables and plots are then generated in real time. For example, each time a new iteration of **simExTargId** takes place updated PCA plots and tables will appear in each sub-directory. The idea is to preserve the reproducibility by recording of parameters used (and session information) and providing an easy to navigate and intuitive sub-directory structure.

The vignette can be viewed here with example raw data acquired on a Thermo FT-ICR mass spectrometer (this requires playing with the creation time of the raw data files to simulate a real-time data collection):

[HTMLvignette](http://bit.ly/2rUQSAk)

**simExTargId** has so far only been tested with Agilent (.d) and Thermo (.RAW/.raw) data files but depending on interest could be readily extended to other instrument manufacturers. 

A technical note describing the package is currently under-review and awaiting decision in the journal *Bioinformatics*.

 
