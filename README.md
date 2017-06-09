# simExTargId (v0.2.0)
Performs simultaneous raw data to mzXML conversion (MSconvert), peak-picking, retention time alignment, grouping (xcms), ESI adducts/in-source fragment artefact identification (CAMERA), peak table output pre-processing, automatic outlier detection by PCA and automatic  statistical analysis dependent on co-variates supplied (MetMSLine), signal drift monitoring and possible MS2 target list identification and real-time shiny-based visualizations are also implemented in real-time during a metabolomic profiling experiment.

The R package utilizes MSConvert (proteoWizard), xcms, CAMERA, MetMSLine, Shiny and many other packages to attempt to implement real-time **sim**ultaneous metabolomic MS1-profiling **Ex**perimentation and statistically relevant MS2 **Targ**et **Id**entification.

A main impetus for development of this package was as an email based early warning system for LC-MS instrumental stoppages/errors but also more subtle changes such as instrument drift and (PCA-based) outlying pooled quality control samples for example. When collecting a dataset of precious/limited samples such as those with a low volume/quantity provided (e.g. mouse sera, dried blood spot punches) it can be particularly poignant if un-noticed instrumental stoppages leads to degradation of your samples or unwanted variation/batch effects.
Instrumental drifts and potentially MS2 targets can be identified in real time using the two shiny apps developed for **simExTargId** *peakMonitor* and *targetId*.

Narrowing the gap to truly "online" target feature MS2 fragmentation is also a primary goal of **simExTargId**. All peak-picking, retention time alignment and grouping is performed by xcms, then ESI adducts and isotopes detected by CAMERA in real time. Following this pre-processing of the peak-table is performed by MetMSLine. Automatic PCA-based outlier detection is performed and a warning email sent to user(s) if a QC sample is outlying. Furthermore, automatic real-time univariate co-variate based statistical analysis is then performed. 
The **simExTargId** workflow can be left running during the MS1-profiling data collection providing a degree of reassurance that serious instrumental difficulties will not go unnoticed and also outlying samples for re-injection and statistically relevant LC-MS feature target can be identified before the end of the experimental worklist.

**SimExTargId** also attempts to inculcate a rigorous organizational directory sub-structure for an experiment and helps to organize data and files for the user in real time. Sub directories, tables and plots are generated in real time. For example each time a new iteration of **simExTargId** takes place

The vignette can be viewed here with example raw data acquired on a Thermo FT-ICR mass spectrometer:

simExTargId has only been tested with Agilent (.d) and Thermo (.RAW/.raw) raw data files but could be relatively easily extended to other instrument manufacturers. 

A technical note describing the package is currently under-review in the journal *Bioinformatics*:

 
