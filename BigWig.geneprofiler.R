#30.04.14

library("ShortRead")
library("rtracklayer")
library("Rsamtools")
library("IRanges")
library("parallel")

f.gene.profiler <- function(f.filename,f.filename.rev="",f.savename,f.annotation,f.sense=TRUE,f.binnumber=30, f.binlength=40, f.upstream= 1000, f.downstream=1000, f.savediskspace=FALSE, f.cores=4) {
  
  # f.filename - BigWig file 
  # f.filename.rev - if empty, f.filename will be used on both strand, if BigWig file, f.filename will be used on fw starnd, f.filename.rev will be used on the reverse strand
  # f.savename - wil be added .profile.Rdat OR .profileS.Rdat OR .profileAS.Rdat
  # f.annotation - Rdat file from "~/Desktop/R/annotation" folder
  # f.sense - if TRUE, sense data will be generated, FALSE -AS data will be genetarted - it is used only if f.filename.rev is not empty
  # f.binnumber - number of bins the feature will be divided
  # f.binlength - length of bins in bp - f.binnumber x f.binnlength = average feature length in bp
  # f.upstream - fragment length collected upstream in bp
  # f.downstream - fragment length collected downstream in bp
  # f.savediskspace=TRUE - will use Rle format in the result file, plotter needs more time to process it, but smaller file size
  # f.cores - no of processor cores
  
  
  # f.annotation: "~/Desktop/R/annotation/S.pombe.EF2/EF2.mRNA.Rdat"
  #               "~/Desktop/R/annotation/S.pombe.2007/sp07.mRNA.Rdat"
  #               "~/Desktop/R/annotation/S.pombe.EF2/EF2.mRNA.nucl.WT040213.Rdat"  OR 260613.Rdat
  #               "~/Desktop/R/annotation/S.pombe.2007/sp07.mRNA.nucl.WT040213.Rdat"  OR 2606613.Rdat
  #               "~/Desktop/R/annotation/N.crassa_or74a_10/Nc10.genes.Rdat"
  #               "~/Desktop/R/annotation/S.pombe.EF2/EF2.introns20to400.Rdat"    
  
###  
  #setwd("~/Documents/seq_data/seq_0008_260613_mononucleosome_Sp_WT_hrp1d_hrp3d_set2d_set1d_set9d_dCD/BigWig_files")
  #f.filename<- "N.sp_st344_WT_polyA_RNAseq_pos.bw"
  #f.filename.rev="N.sp_st344_WT_polyA_RNAseq_neg.bw"
  #f.savename<- "N.sp_st344_WT"
  #f.annotation<-"~/Desktop/R/annotation/S.pombe.EF2/EF2.mRNA.Rdat" 
  #f.sense=TRUE
  #f.binnumber=30
  #f.binlength=40
  #f.upstream= 1000
  #f.downstream=1000
  #f.savediskspace=FALSE
  #f.cores <- 4
###  

  f.excluded=integer(0)
  
  Rle.data <- import(f.filename,as="RleList")               # load data file
  if(f.filename.rev!=""){
    Rle.data.rev<-import(f.filename.rev,as="RleList")       # load reverse data file
    if(f.sense==TRUE){                                # Sense 
      f.save.end<-".profileS.Rdat"
    }else{                                            # Antisense
      f.save.end<-".profileAS.Rdat"
      Rle.data.temp<-Rle.data                         # swapping fw and rev strand information
      Rle.data<-Rle.data.rev
      Rle.data.rev<-Rle.data.temp
      rm(Rle.data.temp)
    }
  }else{
    f.save.end<-".profile.Rdat"
  }
  
  load(f.annotation)                                      # load annotation GRanges object  ann.gr
  
  ######################################## 5' probes
  
  if(f.upstream>0){
    suppressWarnings(temp.gr<-promoters(ann.gr,upstream=f.upstream, downstream=0))    # modified ranges, 5 end, suppress warning about trim
    temp.gr<-trim(temp.gr)                                                            # trim out-of-bound ranges
    
    temp.excluded<-(1:length(temp.gr))[width(temp.gr)==0]                            # if width is 0 (start is 1 larger than end), set start to 1 and width to 1
    start(temp.gr[temp.excluded])<-1
    width(temp.gr[temp.excluded])<-1
    f.excluded<-c(f.excluded,temp.excluded)                                           # add these records to the excluded list
    
    profile5.rle<-RleList(mclapply(1:length(temp.gr), function(x){
      if(as.character(strand(temp.gr[x])=="+")){                                                #fw starnd faetures from fw starnd file
        Rle.data[[as.character(seqnames(temp.gr[x]))]][start(temp.gr[x]):end(temp.gr[x])]
      }else{
        if(f.filename.rev==""){                                                                 #reverse starnd faetures also from the same file if no rev file exist
          Rle.data[[as.character(seqnames(temp.gr[x]))]][end(temp.gr[x]):start(temp.gr[x])]
        }else{
          Rle.data.rev[[as.character(seqnames(temp.gr[x]))]][end(temp.gr[x]):start(temp.gr[x])]     #if reverse file exist, reverse starnd faetures from the rev file 
        }
      }    
    }, mc.cores=f.cores))  
  }else{                                                # if f.upstream=0, we have to define profile.rle
    profile5.rle=RleList()
  }
  
  
  ######################################## feature into f.binnumber bins
  
  profilemid.rle<-RleList(mclapply(1:length(ann.gr), function(x){
    if(as.character(strand(ann.gr[x])=="+")){                                                #fw strand faetures from fw strand file
      f.genestrech(as.numeric(Rle.data[[as.character(seqnames(ann.gr[x]))]][start(ann.gr[x]):end(ann.gr[x])]),f.binnumber) 
    }else{
      if(f.filename.rev==""){                                                                 #reverse strand faetures also from the same file if no rev file exist
        f.genestrech(as.numeric(Rle.data[[as.character(seqnames(ann.gr[x]))]][end(ann.gr[x]):start(ann.gr[x])]),f.binnumber)
      }else{
        f.genestrech(as.numeric(Rle.data.rev[[as.character(seqnames(ann.gr[x]))]][end(ann.gr[x]):start(ann.gr[x])]),f.binnumber)     #if reverse file exist, reverse strand faetures from the rev file 
      }
    }    
  }, mc.cores=f.cores))  
  
  
  ######################################## 3' probes
  
  if(f.downstream>0){
    suppressWarnings(temp.gr<-flank(ann.gr,width=f.downstream, start=FALSE, both=FALSE))    # modified ranges, 3 end suppress warning about trim
    temp.gr<-trim(temp.gr)                                                            # trim out-of-bound ranges
    
    temp.excluded<-(1:length(temp.gr))[width(temp.gr)==0]                            # if width is 0 (start is 1 larger than end), set start to 1 and width to 1
    start(temp.gr[temp.excluded])<-1
    width(temp.gr[temp.excluded])<-1
    f.excluded<-c(f.excluded,temp.excluded)                                           # add these records to the excluded list
    
    profile3.rle<-RleList(mclapply(1:length(temp.gr), function(x){
      if(as.character(strand(temp.gr[x])=="+")){                                                #fw starnd faetures from fw starnd file
        Rle.data[[as.character(seqnames(temp.gr[x]))]][start(temp.gr[x]):end(temp.gr[x])]
      }else{
        if(f.filename.rev==""){                                                                 #reverse starnd faetures also from the same file if no rev file exist
          Rle.data[[as.character(seqnames(temp.gr[x]))]][end(temp.gr[x]):start(temp.gr[x])]
        }else{
          Rle.data.rev[[as.character(seqnames(temp.gr[x]))]][end(temp.gr[x]):start(temp.gr[x])]     #if reverse file exist, reverse starnd faetures from the rev file 
        }
      }    
    }, mc.cores=f.cores))  
  }else{                                                # if f.upstream=0, we have to define profile.rle
    profile3.rle=RleList()
  }
  
  
  ############## combining 5' , feature and 3' rle
  if(f.upstream>0 & f.downstream>0){
    profile.rle<-RleList(mclapply(1:length(profilemid.rle), function(x) c(profile5.rle[[x]],profilemid.rle[[x]],profile3.rle[[x]]),mc.cores=f.cores))
  }else{
    if(f.upstream>0){
      profile.rle<-RleList(mclapply(1:length(profilemid.rle), function(x) c(profile5.rle[[x]],profilemid.rle[[x]]),mc.cores=f.cores))
    }else{
      if(f.downstream>0){
        profile.rle<-RleList(mclapply(1:length(profilemid.rle), function(x) c(profilemid.rle[[x]],profile3.rle[[x]]),mc.cores=f.cores))
      }else{
        profile.rle<-profilemid.rle
      }
    }
  }
  
  
  rle.length<-sum(runLength(profile.rle))             # checking if all features are full length, non full length are excluded and exluded feature numbers are saved in f.excluded
  if(sum(rle.length!=(f.upstream+f.binnumber+f.downstream))>0){
    temp.excluded<-(1:length(ann.gr))[rle.length!=(f.upstream+f.binnumber+f.downstream)]
    f.excluded<-c(f.excluded,temp.excluded)
    f.excluded<-unique(f.excluded)
    print(paste("Lengths anomalie! - Excluded features:", length(f.excluded)))
    print(as.character(mcols(ann.gr[f.excluded])[["Name"]]))
    profile.rle<-profile.rle[-f.excluded]
  }
  

  if(f.savediskspace==FALSE){
    profile.rle<-matrix(unlist(mclapply(1:length(profile.rle),function(x) as.numeric(profile.rle[[x]]), mc.cores=2*f.cores))
                        ,nrow=length(profile.rle), ncol=f.upstream+f.binnumber+f.downstream, byrow=TRUE)  #converting Rle list to matrix of numerical values
  }
  
  save(profile.rle,f.annotation,f.upstream,f.binnumber,f.binlength,f.downstream,f.savediskspace,f.excluded, file=paste(f.savename, f.save.end, sep=""))
  
}


f.genestrech <- function(f.feature,f.binnumber){        # devide an f.feature to f.binnumber bins
  
  #f.binnumber<-as.numeric(f.binnumber)
  f.FElength<- length(f.feature)                        # lenght of the feature
  s.ratio<- f.binnumber/f.FElength                      # theoretical ratio of how many bins will be occupied by one bp of the feature
  a.feature <- rep(s.ratio,f.FElength)                  # s.ratio to every bp of the feature
  a.binnumber <- rep(0.000,f.binnumber)                 # how many bins did we fill up
  FEbins <- rep(0.00, f.binnumber)                      # FEbins will contain the final bins of the feature
  
  
  i <- 1										# ith bp in the feature
  j <- 1										# jth element in the bins
  
  while(a.feature[f.FElength]>0){				# repeat until we didn't divide the feature completely
    if(a.feature[i]>(1-a.binnumber[j])){ 			# if we have to divide the ith bp in the feature 
      FEbins[j]<-FEbins[j]+(1-a.binnumber[j])*f.feature[i]    #FEbins j. elembe az ORF i. elemenek egy resze
      a.feature[i] <-a.feature[i]-(1-a.binnumber[j])              #faeture array i. elemet csokkenteni a resszel amit FEbins-be beletettunk
      if(a.feature[i]<0.0001) a.feature[i]<-0		             #kerekitesi hiba miatt 0 nem mindig 0
      a.binnumber[j]<- 1						             #FEbins j. elem tele van
      j<-j+1								             #ugras az FEbins j+1. elemre - i marad mert a.feature[i] meg >0
    }else{ 									# ha az i. elem osztatlanul megy a j. elembe
      FEbins[j]<-FEbins[j]+(f.feature[i]*a.feature[i])         #FEbins j. elembe az ORF i. elemenek a maradeka
      a.binnumber[j]<- a.binnumber[j]+a.feature[i]               #FEbins array j. elemet novelni a beletett resszel (maximum 1 az eredmeny)
      a.feature[i]<-0							             #ORF i. eleme ures
      i<-i+1								             #ugras az ORF i+1. elemere
      if(a.binnumber[j]==1) j<-j+1			             #ha a FEbins j. eleme tele van, ugras a j+1. elemre
    }
  }
  
  return(FEbins)
}

