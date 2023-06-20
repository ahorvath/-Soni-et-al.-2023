#27.04.15

library("lattice")
library("latticeExtra")
library("caTools")

profile.plotter<-function(f.filename,f.mode="mean",f.quantile=0.5,f.name="NA",f.printmatrix=integer(), f.printmatrix2=integer(), f.xlim=c(0,0), 
                          f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-30,0,30),f.panels=character(0), f.panel.middle=0, f.runmean.segment=TRUE,f.runmean=c(1,1,1),f.oldfilename="",
                          f.finalplot=FALSE, f.log2=FALSE,f.replace0to=0.5,f.convertbackto10=FALSE, f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="",f.mod.annotation="",...){
  
  # f.filename - Rdat filename, created by geneprofiler; if f.filename="" - it will only print the f.printmatrix
  # f.mode="quantile" or "mean",f.quantile=0.5, if f.mode = Quantile, the f.quantile will be generated of every column, otherwise the mean will be calculated (default)
  # if f.log2 = TRUE, data will be converted to log2 scale and geometric averge will be calculated; 0 reads (and < than f.replace0to) will be adjusted to f.replace0to
  # if f.convertbackto10 == TRUE - if f.log2 was TRUE, it willconvert back to decimal before printing (geom average but decimal printing)
  # f.name - name of the series on the graph and in data frame returned from the function
  # f.printmatrix - should be the same length as profile from f.filename, except if adjusted with f.xtrim, new data will be added as last column to f.printmatrix
  #               - if f.filename=="", only f.printmatrix will be printed
  # f.printmatrix2 - if not empty, it should have x values and as many other columns as panels, this will be colored by filling until 0 with grey color, f.xticks should also be given
  #                 - if less colums than panels, only those panels will have this shaded print that has data 
  # f.xtrim(min, max) will trim the data before printing, also returned dataframe will be trimed (x values outside of the range will be lost), if f.xtrim==(0,0), no triming
  # f.xlim(min, max) will set the x limits on the graph, no triming in the data frame
  # f.ylim(min, max) will set the y limits on the graph, no triming in the data frame
  # f.panels = c("Name1", "Name2", ...) - should give name to all columns (including new data), and it will be printed on the panel according the name
  #                                     - if 4 different names, 4 panel will be generated; 1st column(X not included) is the 1st name,  new data is the last name
  #                                     - if more columns than panels, panels will be recycled -e.g. if 2 panels, 1st col is 1st panel, 2nd is 2nd, 3rd col is 1st panel, etc.
  # f.panel.middle - if one number, all middle lines will be at this, if a list, then each panel should have a number where the middle line will be drawn
  # f.xlab; f.ylab  - x and y label
  # f.runmean.segment=TRUE - if TRUE, the segments (such as 5', gene, 3') will be separated and runmean will stop at segment boundaries
  # f.runmean - moving window average, the size of the window in bp - only for printing, the returned data will be the unaveraged, if segments=TRUE, vector with 3 numbers: u1e9pstream, bins, downstream
  # f.oldfilename="" - this has to be provided if f.filename is empty, but we need some of the variables from this file - any of them is ok
  # f.finalplot - TRUE/FALSE - if TRUE, finalpot is given back to the caller, otherwise the f.matrix
  # f.replace0to - used only of f.log2=TRUE, 0 values will be replaced by this number (default is 0.5 - -1 in log2 scale) 
  # f.legendcols - if 0, as many legend columns as the number of columns, otherwise the given number
  # f.middle - only important if no file is loaded, determines the position of the 2nd dotted line, if 0, no line
  # f.mod.annotation - if the Rdat files were generated on the server, you have to add here the same annotation filename/path on your computer!
  # f.xticks - c(a,b,c,d..) - the ticks on x axis 
  # f.filter - a filter file (variable name should be gr.filter) made for the annotation used in the profiler to generate the Rdat file
  
  # ... will be in the par.settings=simpleTheme(); e.g. col=c("green, "red") or other parameters(col, alpha, cex, pch, lty, lwd, font, fill, 
  #                                                                                              border,col.points, col.line, alpha.points, alpha.line)
  # ... e.g. col=c("green","red"); lwd=5
  
  ####
  #f.filter="~/Documents/annotation/S.pombe.EF2/filters/EF2.mRNA.NoHetChr.filter"
  #f.filename<-"exp.hrp1dhrp3d.c043.profileAS.Rdat"
  #f.mode="mean"
  #f.quantile=0.5
  #f.name="NA"
  #f.printmatrix=integer()
  #f.printmatrix2=integer()
  #f.xlim=c(0,0)
  #f.ylim=c(0,0)
  #f.xtrim=c(0,0)
  #f.xticks=c(-30,0,30)
  #f.panels=character(0)
  #f.runmean.segment=TRUE
  #f.runmean=c(1,1,1)
  #f.log2=FALSE
  #f.replace0to=0.5
  #f.xlab="bp"
  #f.ylab="Relative occupancy (log2)"
  #f.legendcols=0
  #f.middle=0
  #f.filter=""
  #f.mod.annotation=""

  ####  
  
  
  if(f.filename !=""){
    load(f.filename)    					# data loaded into profile.rle; annotation file name is in f.annotation
    # profile.rle   - profiles for all features in ann.ger (except f.excluded) in rle format or if f.savediskspace= FALSE in matrix format
    # f.annotation  - name and path of annotation file used in the profile generation; will load ann.gr
    # f.upstream, f.middle, f.downstream - upstream, middle and downstream sequence length in the profile in bp.
    # f.excluded - integer array  of the excluded features in f.annotation file (ann.gr)
    # f.savediskspace - if TRUE, profiles in Rle format, if FALSE, profiles in matrix 
    
    f.middle<-f.binnumber*f.binlength                     # the second dotted line should be f.binnumber*f.binlength
    
    #if(f.mod.annotation==""){
    #load(f.annotation)            # annotation GRanges in ann.gr
    # ann.gr - annotation file used in the profile generation
    #}else{                                                        # if Rdat files were done on an other computer, you have to load the annotation
    #load(f.mod.annotation)
    #}
    
    ################
    # filtering here
    ################
    
    if(f.filter!=""){
      load(f.filter)     # gr.filter is loaded length is the original length of the annotation file
      if(length(f.excluded)>0){
        gr.filter<-gr.filter[-f.excluded]         # excluded values by the profiler should be excluded from the filter 
      }
      profile.rle<-profile.rle[gr.filter,]         # profile.rle is filtered 
    }
    
    print(nrow(profile.rle))
    
    ################
    
    ################
    
    if(f.savediskspace==TRUE){
      profile.num<-matrix(unlist(mclapply(1:length(profile.rle),function(x) as.numeric(profile.rle[[x]]), mc.cores=2*f.cores))
                          ,nrow=length(profile.rle), ncol=f.upstream+f.downstream, byrow=TRUE)  #converting Rle list to matrix of numerical values
      rm(profile.rle)
    }else{                                   # if f.savediskspace = FALSE, profile.rle is already converted to integer matrix
      profile.num<-profile.rle
      rm(profile.rle)
    }
    
    if(f.log2==TRUE){                       # if f.log2 = TRUE, data converted to log2 scale, geometric averge is calculated; 0 reads adjusted to f.replace0to reads
      profile.num[profile.num<f.replace0to]<- f.replace0to     # also smaller than f.replace0to values will be replaced by f.replace0to
      profile.num<-log(profile.num,2)
      profile.num[profile.num==-Inf]<-NA
    }
    
    if(f.mode=="quantile"){
      profile=apply(profile.num,2,function(x) quantile(x,probs=f.quantile, names=FALSE))    #calculate the median of every column (of a geneprofiler data)    
    }else{                         # mean
      profile<-colMeans(profile.num, na.rm=TRUE)
    }
    
    if(f.convertbackto10==TRUE&f.log2==TRUE){   #if we changed it to log2 scale but want to print normal scale (e,g, geom average but print linear scale)
      profile<-2^profile
    }  
  
    if(f.middle==0|f.binlength==1){                        # if no "average gene" feature (only 5 or 3 prime endprofiler) OR binlength=1
      X<-(-f.upstream:(f.middle+f.downstream-1))
    }else{                                  # if bins are used and length !=1
      if(f.upstream>0){
        X<-c(-f.upstream:-1,seq(from=f.binlength/2, to=f.middle-f.binlength/2, by=f.binlength), (f.middle):(f.middle+f.downstream-1))
      }else{
        X<-c(seq(from=f.binlength/2, to=f.middle-f.binlength/2, by=f.binlength), (f.middle):(f.middle+f.downstream-1))
      }
    }
    
    
    avg.data<-as.data.frame(cbind(X,profile))
    colnames(avg.data)<-c("X",f.name)
    if(sum(abs(f.xtrim))>0){               # if f.xtrim()is not (0,0) --> trim the data
      avg.data<-avg.data[avg.data[,1]>=f.xtrim[1]&avg.data[,1]<=f.xtrim[2],]
    } 
  }
  
  ###########################################################################################################################################
  ###############                                if f.printmatrix brings other data series                                    ###############
  ###########################################################################################################################################
  
  if(length(f.printmatrix)>0){                   # if f.printmatrix brings other data series  
    if(sum(abs(f.xtrim))>0){                     # if f.xtrim()is not (0,0) --> trim f.printmatrix
      f.printmatrix<-f.printmatrix[f.printmatrix[,1]>=f.xtrim[1]&f.printmatrix[,1]<=f.xtrim[2],]
    }
    if(f.filename==""){
      avg.data<-f.printmatrix                   # if f.filename=="" - no new data, print only f.printmatrix
    }else{                                                      # otherwise 
      if(nrow(avg.data)!=nrow(f.printmatrix)){                  # check if f.printmatrix length is the same as new profile
        print("X coordinates in new profile don't match X coordintaes in f.printmatrix !!!")
        return(f.printmatrix)
      }else{                                                    # append new profile to f.printmatrix as last column
        while(sum(colnames(f.printmatrix)==f.name)>0){          # if f.printmatrix has already column with f.name, append ".new"
          f.name<-paste(f.name,".new",sep="" )
        }
        avg.data<-cbind(f.printmatrix,avg.data[,2])
        colnames(avg.data)<-c(colnames(f.printmatrix),f.name)
      }
    } 
  }
  # new data plus print.matrix in avg.data dataframe, 1st column is X
  f.printmatrix<-avg.data                                     # f.printmatrix will be handed back to caller
  
  if(f.filename ==""){               # if f.filename is empty, we neeed to load one compatible rle file to get f.upsteam, and other variables
    load(f.oldfilename)    					# data loaded into profile.rle; annotation file name is in f.annotation
  }
  
  if(f.runmean.segment){              # if runmean should stop at segment boundaries
    if(sum(f.runmean>1)>0){
      if(f.upstream>0){
        avg.data[1:f.upstream,]<-as.data.frame(cbind(avg.data[1:f.upstream,1],sapply(2:ncol(avg.data), function(x) runmean(avg.data[1:f.upstream,x],f.runmean[1]))))     #moving window average
        if(length(f.printmatrix2)>0){
          f.printmatrix2[1:f.upstream,2]<- runmean(f.printmatrix2[1:f.upstream,2],f.runmean[1])     #moving window average
        }
      }
      if(f.binnumber>0){
        avg.data[(f.upstream+1):(f.upstream+f.binnumber),]<-as.data.frame(cbind(avg.data[(f.upstream+1):(f.upstream+f.binnumber),1],
                                                                                sapply(2:ncol(avg.data), function(x) runmean(avg.data[(f.upstream+1):(f.upstream+f.binnumber),x],f.runmean[2]))))     #moving window average
        if(length(f.printmatrix2)>0){
          f.printmatrix2[(f.upstream+1):(f.upstream+f.binnumber),2]<- runmean(f.printmatrix2[(f.upstream+1):(f.upstream+f.binnumber),2],f.runmean[2])     #moving window average
        }
      }
      if(f.downstream>0){
        avg.data[(f.upstream+f.binnumber+1):(f.upstream+f.binnumber+f.downstream),]<-as.data.frame(cbind(avg.data[(f.upstream+f.binnumber+1):(f.upstream+f.binnumber+f.downstream),1],
                                                                                sapply(2:ncol(avg.data), function(x) runmean(avg.data[(f.upstream+f.binnumber+1):(f.upstream+f.binnumber+f.downstream),x],f.runmean[3]))))     #moving window average
        if(length(f.printmatrix2)>0){
          f.printmatrix2[(f.upstream+f.binnumber+1):(f.upstream+f.binnumber+f.downstream),2]<- runmean(f.printmatrix2[(f.upstream+f.binnumber+1):(f.upstream+f.binnumber+f.downstream),2],f.runmean[3])     #moving window average
        }
      }
      
      colnames(avg.data)<-colnames(f.printmatrix)
    }
    }else{                            # if runmean should stop go over segment boundaries
      if(f.runmean[1]>1){
        avg.data<-as.data.frame(cbind(avg.data[,1],sapply(2:ncol(avg.data), function(x) runmean(avg.data[,x],f.runmean[1]))))     #moving window average
        colnames(avg.data)<-colnames(f.printmatrix)
        if(length(f.printmatrix2)>0){
          i<-2
          while(i<=length(f.printmatrix2)){
            f.printmatrix2[,i]<- runmean(f.printmatrix2[,i],f.runmean[1])     #moving window average
            i<-i+1
          }
          
        }
      }
  }

 
  ###########################################################################################################################################
  ###############                                            Preparation for printing                                         ###############
  ###########################################################################################################################################
  
  i<-2
  iPanel<-1
  X<-integer(0)
  Reads<-integer(0)
  Panel<-integer(0)
  Experiment<-integer(0)
  
  while(i<=ncol(avg.data)){                                           # converting horizontal data_frame to vertical df 
    X<-c(X,avg.data$X)                                                # Experiment - Panel
    Reads<-c(Reads,avg.data[,i])
    if(length(f.panels)==0){                                           # Panels
      Panel<-c(Panel,rep("Composite plot",nrow(avg.data)))
    }else{
      Panel<-c(Panel,rep(f.panels[iPanel],nrow(avg.data)))
    }
    Experiment<-c(Experiment,rep(colnames(avg.data)[i],nrow(avg.data)))
    i<-i+1
    iPanel<-iPanel+1
    if(iPanel>length(f.panels)){                                   # dataset is distributed between the panels, every new column is to a new panel
      iPanel<-1
    }
  }
  
  print.data<- as.data.frame(cbind(X,Reads,Panel,Experiment))
  print.data$X<-as.numeric(as.character(print.data$X))
  print.data$Reads<-as.numeric(as.character(print.data$Reads))
  
  
  if(f.legendcols==0){            # if f.legendcols not specified, the number of data columns
    f.legendcols<-i-2
  }
  
  ###########################################################################################################################################
  ###############                                            P R I N T I N G                                                  ###############
  ###########################################################################################################################################
  
  if(length(f.printmatrix2)==0){
    final.plot<-xyplot(Reads~X | Panel, data=print.data, groups=Experiment[,drop=TRUE], type="l",par.settings=simpleTheme(...),
                       auto.key=list(columns=f.legendcols, points=FALSE,lines=TRUE),xlab=f.xlab,ylab=f.ylab,
                       scales=list(x=list(relation="same",at=f.xticks),y="same"),
                       panel=function(x,y,...){
                         #panel.xyarea(genemed.matrix2[,1],genemed.matrix2[,2], origin=0,col="lightgrey",col.line="lightgrey",border="lightgrey",lwd=8)
                        panel.xyplot(x,y,...)
                         
                        if(length(f.panel.middle)>1){                                          # if multiple values for f.panel.middle, then use them one by one for the panels
                          if(panel.number()<=length(f.panel.middle)){
                            panel.lines(c(f.panel.middle[panel.number()],f.panel.middle[panel.number()]),c(-100,100), col="black", lwd=2, lty=2)
                          }
                        }else{
                           panel.lines(c(f.panel.middle,f.panel.middle),c(-100,100), col="black", lwd=2, lty=2)    # if only 1 value for f.panel.middle, then use  that for all panels
                        }
                        
                        if(f.middle!=0){
                          panel.lines(c(f.middle,f.middle),c(-100,100), col="black", lwd=2, lty=3)
                        } 
                        
                         
                         panel.lines(c(-10000,10000),c(0,0), col="black", lwd=2, lty=2)        # horizontal line at 0
                       })
  }else{
    final.plot<-xyplot(Reads~X | Panel, data=print.data, groups=Experiment[,drop=TRUE], type="l",par.settings=simpleTheme(...),
                       auto.key=list(columns=f.legendcols, points=FALSE,lines=TRUE),xlab=f.xlab,ylab=f.ylab,
                       scales=list(x=list(relation="same",at=f.xticks),y="same"),
                       
                       panel=function(x,y,...){
                
                         if(panel.number()<=(ncol(f.printmatrix2)-1)){
                           panel.xyarea(f.printmatrix2[,1],f.printmatrix2[,(panel.number()+1)], origin=0,col="lightgrey",col.line="lightgrey",border="lightgrey",lwd=8)
                         }
                         
                         panel.xyplot(x,y,...)
                         
                         if(length(f.panel.middle)>1){                                          # if multiple values for f.panel.middle, then use them one by one for the panels
                           if(panel.number()<=length(f.panel.middle)){
                             panel.lines(c(f.panel.middle[panel.number()],f.panel.middle[panel.number()]),c(-100,100), col="black", lwd=2, lty=2)
                           }
                         }else{
                           panel.lines(c(f.panel.middle,f.panel.middle),c(-100,100), col="black", lwd=2, lty=2)    # if only 1 value for f.panel.middle, then use  that for all panels
                         }
                         
                         if(length(f.panel.middle)>1){                                          # if multiple values for f.panel.middle, then use them one by one for the panels
                           if(panel.number()<=length(f.panel.middle)){
                             panel.lines(c(f.panel.middle[panel.number()],f.panel.middle[panel.number()]),c(-100,100), col="black", lwd=2, lty=2)
                           }
                         }else{
                           panel.lines(c(f.panel.middle,f.panel.middle),c(-100,100), col="black", lwd=2, lty=2)    # if only 1 value for f.panel.middle, then use  that for all panels
                         }
                         
                         if(f.middle!=0){
                           panel.lines(c(f.middle,f.middle),c(-100,100), col="black", lwd=2, lty=3)
                         } 
                         
                       })
  }
  
  
  
  if(sum(abs(f.ylim))>0){                               # if ylim was specified
    final.plot<-update(final.plot, ylim=f.ylim, evaluate=FALSE)
  }
  if(sum(abs(f.xlim))>0){                               # if xlim was specified
    if(sum(abs(f.ylim))>0){                               # if ylim was specified
      final.plot<-update(final.plot, ylim=f.ylim, xlim=f.xlim, evaluate=FALSE)
    }else{
      final.plot<-update(final.plot, xlim=f.xlim, evaluate=FALSE)
    }
  }
  
  #final.plot<-update(final.plot, panel.xyarea(f.printmatrix2[,1],f.printmatrix2[,2], origin=0,col="lightgrey",col.line="lightgrey",border="lightgrey",lwd=8))
  #final.plot<-final.plot+ layer(panel.xyarea(print.data.2$X.2, print.data.2$Reads.2, groups=print.data.2$Panel.2, origin=0,col="lightgrey",col.line="lightgrey",border="lightgrey",lwd=8))
  #final.plot<-final.plot+ layer(panel.xyarea(f.printmatrix2[,1],f.printmatrix2[,2], origin=0,col="lightgrey",col.line="lightgrey",border="lightgrey",lwd=8))
  
  print(final.plot)
  
  # update(final.plot,
  #       panel = function(...,){
  #        panel.xyplot(...)
  #       panel.lines(...)
  #    })
  #final.plot<-xyplot(Reads~X | Panel, data=print.data, groups=Experiment[,drop=TRUE],type=c("p","p","p","l"),
  #                   distribute.type=TRUE, ylim=c(ylim.min,ylim.max),col=c("black","black","black","red"),
  #                  lty=c(1,1,1,1), pch=c(16,45,45,1),cex=c(1,2,2,1),lwd=c(1,1,1,2), as.table=TRUE,
  #                 scales=list(x=list(relation="same",at=c(0,6,12,18,24,30,36,42,48)),y="same"),
  #                xlab="time (h)", ylab="Expression Ratio (log2)",
  #               main= paste(f.name,"|",g.name,"| Phase:", mean.sinphase,"| SinCor2:", mean.SinCor2,"| SinCor3:",mean.SinCor3,sep=""))
  
  #trellis.focus("panel", column=2, row=1)
  #panel.lines(c(0,0),c(-100,100), col="gray", lwd=2, lty=2)                               # dotted line at 0
  #trellis.unfocus()
  
  
  ###########################################################################################################################################
  ###############                                              R E T U R N                                                    ###############
  ########################################################################################################################################### 
  
  
  if(f.finalplot==TRUE){
    return(final.plot)
  }else{
    return(f.printmatrix)                               # f.printmatrix has previous f.printmatrix with new trim and new profile in the last column 
  }
}  

#X[peaks(profile,span=99)]

peaks<-function(series,span=3, ties.method = "first")
{
  if((span <- as.integer(span)) %% 2 != 1) stop("'span' must be odd")
  z <- embed(series, span)
  s <- span%/%2
  v <- max.col(z, ties.method=ties.method) == 1 + s
  pad <- rep(FALSE, s)
  result <- c(pad, v, pad)
  result
}