#' xcms (find by formula) results viewer shiny application
#' @param xcmsOutput character full path to xcms results file if argument is not supplied a tcltk window will open
#' @param launchApp logical default is true. If true peakMonitor shiny app will open in your default web browser.
#' @export
peakMonitor <- function(xcmsOutput=NULL, obsNames=NULL, launchApp=TRUE){
  # error handling
  if(is.null(obsNames)){
    stop('The names of the observation columns must be supplied')
  }
  
  if(is.null(xcmsOutput)){
    xcmsOutput <- tcltk::tclvalue(tcltk::tkgetOpenFile(title="select your xcms results file"))
  }
  
  if(is.character(xcmsOutput)){
  if(!grepl('\\.csv$', xcmsOutput)){
    stop('xcmsOutput file must be a .csv')
  }
    xcmsOutput <- as.data.frame(data.table::fread(xcmsOutput, header=TRUE, stringsAsFactors = FALSE))
  }
  
  if(!is.data.frame(xcmsOutput)){
    stop("xcmsOutput must be supplied as either a full path to an xcms output file or as a data.frame")
  }
  # GLOBAL
  
  
  skiptoptwo_rows<-all_content[-c(1,2,length(all_content))]
  xcms_results<-read.csv(textConnection(skiptoptwo_rows), header = TRUE,row.names=NULL, stringsAsFactors = FALSE)
  xcms_results<-xcms_results[,-1]
  xcms_colnames<-colnames(xcms_results)[2:ncol(xcms_results)]
  xcms_results<-xcms_results[,-ncol(xcms_results)]
  colnames(xcms_results)<-xcms_colnames
  ###remove columns containing only NAs
  
  xcms_results<-xcms_results[,colSums(is.na(xcms_results)) != nrow(xcms_results)]
  ###replace NAs remaining with zeros
  xcms_results[is.na(xcms_results)]<-0
  NamesColumn<-colnames(xcms_results)[grep("Name|Notes",colnames(xcms_results))][1]
  xcms_results[,NamesColumn]<-gsub("-",".",xcms_results[,NamesColumn])
  
  xcms_results<-xcms_results[order(xcms_results$File),]
  unique.compound.name<-unique(xcms_results[,NamesColumn])
  unique.sample.names<-unique(xcms_results$File)
  
  xcms_results_colnames<-c("Height","Area","Abund","RT","Mass..Tgt.",colnames(xcms_results)[grep("Diff..ppm.|Diff..Tgt..ppm.",colnames(xcms_results))])
  ###reshape dataframe so that formula appear in columns
  xcms_results_reshaped<-lapply(c(1:length(unique.compound.name)), function (x)
    xcms_results[which(xcms_results[,NamesColumn] == unique.compound.name[x]),c(xcms_results_colnames,"File","Score")])
  ###sort by height or score? prior to removal of duplicate file names
  xcms_results_reshaped<-lapply(xcms_results_reshaped,function(x) {x<-x[order(x$Height,decreasing=T),]
                                                                 x<-x[duplicated(x$File)==F,]
                                                                 x<-x[order(x$File),]})
  
  xcms_results_reshaped<-lapply(c(1:length(xcms_results_reshaped)),function(x) {
    lapply(c(1:length(unique.sample.names)),function(y) {
      if((unique.sample.names %in% xcms_results_reshaped[[x]][,"File"])[[y]]==T){
        sample.result<-xcms_results_reshaped[[x]][which(xcms_results_reshaped[[x]]$File %in% unique.sample.names[y]),]
        return(sample.result)
      } else {
        dummy.m<-data.frame(matrix(0,ncol=ncol(xcms_results_reshaped[[1]]),nrow=1))
        colnames(dummy.m)<-colnames(xcms_results_reshaped[[1]])
        dummy.m$File<-unique.sample.names[y]
        return(dummy.m)
      }
    })
  })
  
  xcms_results_reshaped<-lapply(xcms_results_reshaped,function(x) do.call(rbind,x))
  
  xcms_results_reshaped<-do.call(cbind,xcms_results_reshaped)
  
  xcms_results_reshaped<-xcms_results_reshaped[,-grep("File|Score",colnames(xcms_results_reshaped))]
  ###bind the sample names on to the reshaped dataframe
  row.names(xcms_results_reshaped)<-unique.sample.names
  ###add column names####
  colnames(xcms_results_reshaped)<-paste(rep(unique.compound.name,each=length(xcms_results_colnames)),"_",xcms_results_colnames,sep="")
  ###identify run order number from file name and regular expression
  runorder<-regexpr("_[0123456789]*_",row.names(xcms_results_reshaped))
  string.length<-attr(runorder,"match.length")
  ###for loop to extract run order information from sample names
  for (i in 1:nrow(xcms_results_reshaped))
  {
    unique.sample.names[i]<-gsub("_","",substr(unique.sample.names[i],as.numeric(runorder[i]),as.numeric(runorder[i]+string.length[i]-1)))
  }
  
  xcms_results_reshaped<-xcms_results_reshaped[order(as.numeric(unique.sample.names)),]
  
  ###function to create sortable numbers in file names
  pad_int<-function(n,scale){
    out_string<-paste(10*scale + n,sep='')
    out_string<-substr(out_string,2,nchar(out_string))
    return(out_string)
  }
  
  ###loop to replace sample number and injection order number with padded numbers
  for(i in 1:nrow(xcms_results_reshaped)){
    ###find number locations
    number.locations<-gregexpr("_[0123456789]*_",row.names(xcms_results_reshaped)[i])
    if(number.locations[[1]][1]!=(-1))
    {
      number.locations<-c(number.locations,gregexpr("_[0123456789]*.d$",row.names(xcms_results_reshaped)[i]))
      
      number.strings<-sapply(c(1:length(number.locations)),function(x) substr(row.names(xcms_results_reshaped)[i],number.locations[[x]],(as.numeric(number.locations[[x]])+attr(number.locations[[x]],"match.length")-1)))
      row.names(xcms_results_reshaped)[i]<-gsub(number.strings[1],paste("_",pad_int(as.numeric(gsub("_","",number.strings[1])),1000),"_",sep=""),row.names(xcms_results_reshaped)[i])                                                      
      if(length(number.strings)>1){
        row.names(xcms_results_reshaped)[i]<-gsub(number.strings[2],paste("_",pad_int(as.numeric(gsub(".d","",gsub("_","",number.strings[2]))),1000),".d",sep=""),row.names(xcms_results_reshaped)[i])                                                      
      }
    }
  }
  
  xcms_results_colnames<-xcms_results_colnames[-grep("Diff..ppm.|Diff..Tgt..ppm.|Mass..Tgt.",xcms_results_colnames)]
  
  Log2.zero.df<-xcms_results_reshaped[,grep("Area",colnames(xcms_results_reshaped),ignore.case=T)]
  
  Log2.zero.df[Log2.zero.df==0]<-min(unlist(Log2.zero.df)[unlist(Log2.zero.df)!=0])
  ###calculate PCA plot of raw data and plot using rCharts
  PCA.result<-pca(apply(Log2.zero.df,2,log2), nPcs = 2, method="svd",scale="pareto", centre=TRUE, cv="q2")
  
  ##Get the estimated scores
  scores<-as.data.frame(PCA.result@scores)
  
  ###generate rCharts plot of data
  scores$Sample_type<-grepl("QC",row.names(scores))
  scores$Sample_type<-ifelse(scores$Sample_type==TRUE,"QC","Sample")
  scores$Sample_type[grep("blank",row.names(scores),ignore.case=T)]<-"Blank"
  scores$Sample_type[grep("Global",row.names(scores),ignore.case=T)]<-"Global_QC"
  ###add sample name column
  scores$Sample_name<-row.names(scores)
  
  ##Get the loadings
  loadings<-as.data.frame(PCA.result@loadings)
  loadings$Metabolite_PeakArea_name<-row.names(loadings)
  ###calculate mean all compounds
  VariableTypes<-unique(gsub(".+_","",colnames(xcms_results_reshaped)))
  MeanAllCompounds.df<-sapply(VariableTypes,function(x){ apply(xcms_results_reshaped[,grep(x,colnames(xcms_results_reshaped))],1,mean)})
  colnames(MeanAllCompounds.df)<-paste0("MeanAllCompounds_",VariableTypes)
  MeanAllCompounds.df[,"MeanAllCompounds_Mass..Tgt."]<-mean(as.numeric(MeanAllCompounds.df[,"MeanAllCompounds_Mass..Tgt."]))
  xcms_results_reshaped<-cbind(xcms_results_reshaped,MeanAllCompounds.df)
  # SERVER
  xcmsserver <- shinyServer(function(input, output){ 
    
    output$PCA_plot_1 <- renderChart({
      ###print model summary
      #print(xtable(data.frame(PCA.result@R2,PCA.result@R2cum,PCA.result@cvstat)),type="html")  
      m1<-rPlot(x="PC1",y="PC2",color="Sample_type",type="point",data=scores,size = list(const = 5))
      
      m1$set(dom="PCA_plot_1")
      
      return(m1)
      
    })
    
    output$PCA_plot_2 <- renderChart({
      ###print model summary
      #print(xtable(data.frame(PCA.result@R2,PCA.result@R2cum,PCA.result@cvstat)),type="html")  
      m2<-rPlot(x="PC1",y="PC2",color="Sample_name",type="point",data=scores,size = list(const = 5))
      
      m2$set(dom="PCA_plot_2")
      
      #m2$params$width<-600
      #m2$params$height<-500
      
      return(m2)
      
    })
    
    output$PCA_plot_loadings <- renderChart({
      ###print model summary
      #print(xtable(data.frame(PCA.result@R2,PCA.result@R2cum,PCA.result@cvstat)),type="html")  
      m3<-rPlot(x="PC1",y="PC2",color="Metabolite_PeakArea_name",type="point",data=loadings,size = list(const = 5))
      
      m3$set(dom="PCA_plot_loadings")
      
      #m2$params$width<-600
      #m2$params$height<-500
      
      return(m3)
      
    })
    
    output$xcms_results_plot <- renderChart({
      
      if(!is.null(input$Exclude))
      {
        plot.df<-xcms_results_reshaped[-which(row.names(xcms_results_reshaped) %in% input$Exclude),grep(input$Compound,colnames(xcms_results_reshaped))]
        
      }else{
        ###subset by compound
        plot.df<-xcms_results_reshaped[,grep(input$Compound,colnames(xcms_results_reshaped))]
      }
      ###subset by result type
      plot.df<-data.frame(cbind(row.names(plot.df),plot.df[,grep(input$Result_type,colnames(plot.df))]))
      colnames(plot.df)<-c("Sample",input$Result_type)
      ###calculate the percentage of the mean
      if(input$Deviation_mean=="YES")
      {
        plot.df[,2]<-(as.numeric(as.character(plot.df[,2]))/mean(as.numeric(as.character(plot.df[,2]))))*100 
      }
      
      if(length(grep("QC",plot.df[,1]))>0)
      {
        plot.df$avg <- mean(as.numeric(as.character(plot.df[grep("QC",plot.df[,1]),2])))#,length.out=nrow(plot.df))
        plot.df$sdplus <- plot.df$avg + sd(as.numeric(as.character(plot.df[grep("QC",plot.df[,1]),2])))
        plot.df$sdneg <-  plot.df$avg - sd(as.numeric(as.character(plot.df[grep("QC",plot.df[,1]),2])))
      } else {
        plot.df$avg <- mean(as.numeric(as.character(plot.df[,2])))#,length.out=nrow(plot.df))
        plot.df$sdplus <- plot.df$avg + sd(as.numeric(as.character(plot.df[,2])))
        plot.df$sdneg <-  plot.df$avg - sd(as.numeric(as.character(plot.df[,2])))
        
      }
      
      ###add colour for QC vs sample
      plot.df$Sample_type<-grepl("QC",plot.df[,1])
      plot.df$Sample_type<-ifelse(plot.df$Sample_type==TRUE,"QC","Sample")
      plot.df$Sample_type[grep("blank",plot.df[,1],ignore.case=T)]<-"Blank"
      plot.df$Sample_type[grep("Global",plot.df[,1],ignore.case=T)]<-"Global_QC"
      
      m4<-rPlot(x="Sample",y=colnames(plot.df)[2],color="Sample_type",type="bar",data=plot.df,size = list(const = 5))
      
      #m1$set(pointSize = 0, lineWidth = 1)
      
      m4$set(dom="xcms_results_plot")
      m4$layer(y='avg', copy_layer=T, type='line', color=list(const='red'))
      m4$layer(y='sdplus', copy_layer=T, type='line', color=list(const='green'))
      m4$layer(y='sdneg', copy_layer=T, type='line', color=list(const='green'))
      
      return(m4)
      
    })
    ###return the coefficient of variance
    output$QC_CV_text<-renderText({
      
      if(!is.null(input$Exclude))
      {
        plot.df<-xcms_results_reshaped[-which(row.names(xcms_results_reshaped) %in% input$Exclude),grep(input$Compound,colnames(xcms_results_reshaped))]
        
      }else{
        ###subset by compound
        plot.df<-xcms_results_reshaped[,grep(input$Compound,colnames(xcms_results_reshaped))]
      }
      plot.df<-data.frame(cbind(row.names(plot.df),plot.df[,grep(input$Result_type,colnames(plot.df))]))
      
      if(length(grep("QC",plot.df[,1]))>0)
      {
        ###only include QC samples into calculation
        QC.df<-plot.df[grep("QC",plot.df[,1]),]
        
        if(length(grep("Global",QC.df[,1],ignore.case=T))!=nrow(QC.df))
        {
          if(length(grep("Global",QC.df[,1],ignore.case=T))>0)
          {
            ##if not only global QCs then remove global QCs from calc
            QC.df<-QC.df[-grep("Global",QC.df[,1],ignore.case=T),] 
          }
        }
        
        ###calculate QC average and CV%
        QC.avg <- mean(as.numeric(as.character(QC.df[,2])))#,length.out=nrow(plot.df))
        QC.cv <- (sd(as.numeric(as.character(QC.df[,2])))/QC.avg)*100
        
        QC_CV_text<-paste("The current CV% for the ",input$Result_type," of ",input$Compound," in the quality controls is ",prettyNum(round(QC.cv,digits=3),big.mark=",",scientific=F),sep="")
        
        return(QC_CV_text)
      } else {
        
        QC_CV_text<-paste0("No Quality control samples were identified to calculate the QC CV% for the ",input$Result_type)
        return(QC_CV_text) 
      }
      
    })
    
    output$Sample_CV_text<-renderText({
      if(!is.null(input$Exclude))
      {
        plot.df<-xcms_results_reshaped[-which(row.names(xcms_results_reshaped) %in% input$Exclude),grep(input$Compound,colnames(xcms_results_reshaped))]
        
      }else{
        ###subset by compound
        plot.df<-xcms_results_reshaped[,grep(input$Compound,colnames(xcms_results_reshaped))]
      }
      plot.df<-data.frame(cbind(row.names(plot.df),plot.df[,grep(input$Result_type,colnames(plot.df))]))
      
      if(length(which(grepl("QC",plot.df[,1])==F))==0)
      {
        return(paste0("No non-QC samples were identified to calculate the sample CV% for the ",input$Result_type))
        
      } else {
        ###only include Samples into calculation
        Sample.df<-plot.df[which(grepl("QC",plot.df[,1])==F),]
        
        ###remove blanks
        if(length(grep("blank",Sample.df[,1],ignore.case=T))==nrow(Sample.df))
        { 
          return(paste0("No non-QC samples were identified to calculate the sample CV% for the ",input$Result_type))
          
        } else {
          
          if(length(grep("blank",Sample.df[,1],ignore.case=T))!=0)
          {
            ##if not only blanks then remove blanks from calc
            Sample.df<-Sample.df[-grep("blank",Sample.df[,1],ignore.case=T),]
          } 
          
          if(sum(as.numeric(as.character(Sample.df[,2])))==0)
          {
            
            Sample_CV_text<-paste0("No non-zero value samples were detected for the ",input$Result_type)
            return(Sample_CV_text) 
            
          } else {
            ###calculate Sample average and CV%
            Sample.avg <- mean(as.numeric(as.character(Sample.df[,2])))#,length.out=nrow(plot.df))
            Sample.cv <- (sd(as.numeric(as.character(Sample.df[,2])))/Sample.avg)*100
            
            Sample_CV_text<-paste("The current CV% for the ",input$Result_type," of ",input$Compound," in the samples is ",prettyNum(round(Sample.cv,digits=3),big.mark=",",scientific=F),sep="")
            
            return(Sample_CV_text)
            
          }
        }
      }
    })
    
    
    
    ###return the average
    output$Av_text<-renderText({
      if(!is.null(input$Exclude))
      {
        plot.df<-xcms_results_reshaped[-which(row.names(xcms_results_reshaped) %in% input$Exclude),grep(input$Compound,colnames(xcms_results_reshaped))]
        
      }else{
        ###subset by compound
        plot.df<-xcms_results_reshaped[,grep(input$Compound,colnames(xcms_results_reshaped))]
      }
      plot.df<-data.frame(cbind(row.names(plot.df),plot.df[,grep(input$Result_type,colnames(plot.df))]))
      
      if(length(grep("QC",plot.df[,1]))>0)
      {
        ###only include QC samples into calculation
        QC.df<-plot.df[grep("QC",plot.df[,1]),]
        if(length(grep("Global",QC.df[,1],ignore.case=T))!=nrow(QC.df))
        {
          if(length(grep("Global",QC.df[,1],ignore.case=T))>0)
          {
            ##if not only global QCs then remove global QCs from calc
            QC.df<-QC.df[-grep("Global",QC.df[,1],ignore.case=T),] 
          }
        }
        
        avg <- mean(as.numeric(as.character(QC.df[,2])))#,length.out=nrow(plot.df))
        
        Av_text<-paste("The current average ",input$Result_type," of ",input$Compound," in the quality controls is ",prettyNum(round(avg,digits=3),big.mark=",",scientific=F),sep="")
        
        return(Av_text)
      } else {
        Av_text<-paste0("No Quality control samples were identified to calculate the average QC ",input$Result_type)
        return(Av_text)  
      }
    })
    
    output$Sample_Av_text<-renderText({
      
      if(!is.null(input$Exclude))
      {
        plot.df<-xcms_results_reshaped[-which(row.names(xcms_results_reshaped) %in% input$Exclude),grep(input$Compound,colnames(xcms_results_reshaped))]
        
      }else{
        ###subset by compound
        plot.df<-xcms_results_reshaped[,grep(input$Compound,colnames(xcms_results_reshaped))]
      }
      plot.df<-data.frame(cbind(row.names(plot.df),plot.df[,grep(input$Result_type,colnames(plot.df))]))
      
      if(length(which(grepl("QC",plot.df[,1])==F))>0)
      {
        Sample.df<-plot.df[which(grepl("QC",plot.df[,1])==F),]
        ###remove blanks
        if(length(grep("blank",Sample.df[,1],ignore.case=T))!=nrow(Sample.df))
        {
          if(length(grep("blank",Sample.df[,1],ignore.case=T))!=0)
          {  
            ##if not only blanks then remove blanks from calc
            Sample.df<-Sample.df[-grep("blank",Sample.df[,1],ignore.case=T),]
          }
        }
        
        ###calculate Sample average and CV%
        Sample.avg <- mean(as.numeric(as.character(Sample.df[,2])))#,length.out=nrow(plot.df))
        
        Sample_Av_text<-paste("The current average ",input$Result_type," of ",input$Compound," in the samples is ",prettyNum(round(Sample.avg,digits=3),big.mark=",",scientific=F),sep="")
        return(Sample_Av_text) 
      } else {
        return(paste0("No non-QC samples were identified to calculate the average sample ",input$Result_type))
      }
    })
    
    output$Mass_accuracy<-renderText({
      
      if(!is.null(input$Exclude))
      {
        plot.df<-xcms_results_reshaped[-which(row.names(xcms_results_reshaped) %in% input$Exclude),grep(input$Compound,colnames(xcms_results_reshaped))]
        
      }else{
        ###subset by compound
        plot.df<-xcms_results_reshaped[,grep(input$Compound,colnames(xcms_results_reshaped))]
      }
      plot.df<-data.frame(cbind(row.names(plot.df),plot.df[,grep("Diff..ppm.|Diff..Tgt..ppm.|Mass..Tgt.",colnames(plot.df))]))
      
      tgtmass<-plot.df[,grep("Mass..Tgt.",colnames(plot.df))]
      
      tgtmass<-unique(tgtmass[tgtmass!=0])
      
      ppm <- mean(abs(as.numeric(as.character(plot.df[,grep("Diff..ppm.|Diff..Tgt..ppm.",colnames(plot.df))]))))#,length.out=nrow(plot.df))
      
      ppmsd <- sd(abs(as.numeric(as.character(plot.df[,grep("Diff..ppm.|Diff..Tgt..ppm.",colnames(plot.df))]))))
      
      ppm_text<-paste("The current average ppm mass accuracy of ",input$Compound," (expected mass: ",tgtmass,") in all samples is ",prettyNum(round(ppm,digits=3),big.mark=",",scientific=F)," (+/- ",prettyNum(round(ppmsd,digits=3),big.mark=",",scientific=F)," s.d.)",sep="")
      
      return(ppm_text)
    })
    
    output$Retention_time<-renderText({
      if(!is.null(input$Exclude))
      {
        plot.df<-xcms_results_reshaped[-which(row.names(xcms_results_reshaped) %in% input$Exclude),grep(input$Compound,colnames(xcms_results_reshaped))]
        
      }else{
        ###subset by compound
        plot.df<-xcms_results_reshaped[,grep(input$Compound,colnames(xcms_results_reshaped))]
      }
      plot.df<-data.frame(cbind(row.names(plot.df),plot.df[,grep("RT",colnames(plot.df))]))
      
      plot.df<-plot.df[plot.df[,2]!=0,]
      
      RTavg<-mean(as.numeric(as.character(plot.df[,2])))#,length.out=nrow(plot.df))
      
      RTsd<-sd(as.numeric(as.character(plot.df[,2])))
      
      RT_text<-paste("The current average Retention time for ",input$Compound," in all samples is ",prettyNum(round(RTavg,digits=3),big.mark=",",scientific=F)," (+/- ",prettyNum(round(RTsd,digits=3),big.mark=",",scientific=F)," s.d.) ",sep="")
      
      return(RT_text)
    })
    
    
  })
  
  # UI
  xcmsui <- shinyUI(fluidPage(
    titlePanel(paste0("xcms results app (",basename(xcmsOutput),")")),
    sidebarLayout(
      
      sidebarPanel(h4("xcms results function options:"),
                   
                   selectInput("Compound",label=tags$b("Select the name of the compound to display :"),choices=c("MeanAllCompounds",sort(unique.compound.name))),
                   selectInput("Result_type",label=tags$b("Select the xcms result type to display :"),choices=c(xcms_results_colnames)),
                   selectInput("Exclude",label=tags$b("Select any samples you wish to exclude:"),choices=c(row.names(xcms_results_reshaped)),multiple=T),
                   radioButtons("Deviation_mean",label=tags$b("Percentage of the mean:"),choices=c("NO","YES"))              
      ),
      
      mainPanel(h4("xcms results plot"),width=8,
                showOutput("xcms_results_plot","polycharts"),
                h5(strong("CV%")),
                (textOutput("QC_CV_text")),
                (textOutput("Sample_CV_text")),
                h5(strong("Average")),
                (textOutput("Av_text")),
                (textOutput("Sample_Av_text")),
                h5(strong("Mass accuracy")),
                (textOutput("Mass_accuracy")),
                h5(strong("Retention time")),
                (textOutput("Retention_time")),
                tags$br(),
                h5(strong("PCA all features (scores)")),
                showOutput("PCA_plot_1","polycharts"),
                showOutput("PCA_plot_2","polycharts"),
                h5(strong("PCA all features (loadings)")),
                showOutput("PCA_plot_loadings","polycharts")
      )
      
    )
    
  ))
  
  shiny::runApp(list(ui = xcmsui, server = xcmsserver), 
                launch.browser = browserLaunch)
  
} # end function
