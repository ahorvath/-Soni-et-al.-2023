# GeneProfiler
library(GenomicRanges)
library(rtracklayer)
library(limma)
library(ggplot2)

fw.bw.files <- dir(pattern="*.f.bw$", path = "bigwigs", full.names = T)
rev.bw.files <- dir(pattern="*.r.bw$", path = "bigwigs", full.names = T)
track.names <- gsub("...bw$", "", apply(strsplit2(basename(fw.bw.files), split = "\\.")[,1,drop = F], 1, paste0, collapse = "_"))
f.annotation="annotation/S.pombe.EF2/EF2.exons25to5123.Rdat"
load(f.annotation)

source("BigWig.geneprofiler.R")
source("profile.plotter.R")
source("infer_average_profile.R")
binnumber <- 40
lapply(1:length(fw.bw.files), function(x) f.gene.profiler(fw.bw.files[[x]], rev.bw.files[[x]], track.names[[x]],f.annotation,f.sense=TRUE,f.binnumber=binnumber, f.binlength=1,f.upstream=0, f.downstream=0,f.cores=14))
#lapply(1:length(fw.bw.file.list1), function(x) f.gene.profiler(fw.bw.file.list1[[x]], rev.bw.file.list2[[x]], file.list3[[x]],f.annotation,f.sense=FALSE,f.binnumber=binnumber, f.binlength=1,f.upstream=0, f.downstream=0,f.cores=14))

S.profiles <- dir(pattern="^sp.*profileS.Rdat$") # exon profile
AS.profiles <- dir(pattern="^sp.*profileAS.Rdat$") # exon profile

source("BigWig.normalize.R")
lapply(1:length(fw.bw.files), function(x) f.BW.norm(fw.bw.files[[x]], f.filename2=rev.bw.files[[x]],f.mode="quantile",f.normvalue=175,f.chr="",f.normfile=S.profiles[[x]], f.quantile=0.5))
#lapply(1:length(fw.bw.files), function(x) f.BW.norm(fw.bw.files[[x]], f.filename2=rev.bw.files[[x]],f.mode="quantile",f.normvalue=200,f.chr="",f.normfile=AS.profiles[[x]], f.quantile=0.5))

norm.fw.bw.file.list1<-dir(pattern="^N.*.f.bw$")
norm.rev.bw.file.list2<-dir(pattern="^N.*.r.bw$")
# N bigwig files
norm.file.list3<-gsub("...bw$","", basename(norm.fw.bw.file.list1))
#f.annotation="annotation/aggregated_introns_AH_Chen_et_al_2018.Rdat"
f.annotation <- "aggregated_introns_AH_Burke_et_al_2018.Rdat"
lapply(1:length(norm.fw.bw.file.list1), function(x) f.gene.profiler(norm.fw.bw.file.list1[[x]], norm.rev.bw.file.list2[[x]], norm.file.list3[[x]],f.annotation,f.sense=TRUE,f.binnumber=binnumber, f.binlength=1,f.upstream=20, f.downstream=20,f.cores=14))

norm.S.profiles <- dir(pattern="^N.*profileS.Rdat$") # exon profile
#norm.AS.profiles <- dir(pattern="^N.*profileAS.Rdat$") # exon profile

# Profile plotter
#ser.order <- c(3,2,1)
ordered.norm.S.profiles <- norm.S.profiles[ser.order]
ser.name <- c("WT_Rep1", "ctr1d_Rep2", "ctr1d_Rep1")
cols <- c( "lightblue", "darkblue", "grey")
print.S <- integer()
normalised.profiles <-dir(pattern = "^N.*profileS.Rdat$") # intron profile
for (i in 1:length(ser.order)) {
  print(cols[i])
  print.S<-profile.plotter(ordered.norm.S.profiles[i],f.mode="mean",f.quantile=0.5,f.name=ser.name[i],f.printmatrix=print.S, f.printmatrix2=integer(), f.xlim=c(0,0),
                           f.ylim=c(0,0),f.xtrim=c(0,0),f.xticks=c(-20, 0, 40, 60),f.panels=character(0), f.runmean.segment=TRUE,f.runmean=c(1,1,1),f.log2=T,f.replace0to=1,
                           f.xlab="bp", f.ylab="Relative occupancy (log2)", f.legendcols=0, f.middle=0,f.filter="",f.mod.annotation="", col=cols,lwd=c(3,3,3), lty=c(1,1,1))
}

source("infer_average_profile.R")

mean.table <- infer.average.profile(ordered.norm.S.profiles, sample.names = ser.name, f.fun ="median", f.binnumber = 40, f.upstream = 20, f.downstream = 20)
load(f.annotation)
ann.median.table <- cbind(data.frame(ann.gr)[,c(1:3,6,4,5,8,7,9:36)], mean.table[,c(1:9)])
ann.median.table <- transform(ann.median.table, Intron_relative_abundance_WT_Rep1 = WT_Rep1.body/((WT_Rep1.upstream + WT_Rep1.downstream)/2),
                                                Intron_relative_abundance_ctr1d_Rep1 = ctr1d_Rep1.body/((ctr1d_Rep1.upstream + ctr1d_Rep1.downstream)/2),
                                                Intron_relative_abundance_ctr1d_Rep2 = ctr1d_Rep2.body/((ctr1d_Rep2.upstream + ctr1d_Rep2.downstream)/2)
                              )
ann.median.table <- transform(ann.median.table, Intron_relative_abundance_ctr1d_Rep1_minus_WT_Rep1 = Intron_relative_abundance_ctr1d_Rep1 - Intron_relative_abundance_WT_Rep1,
                                                Intron_relative_abundance_ctr1d_Rep2_minus_WT_Rep1 = Intron_relative_abundance_ctr1d_Rep2 - Intron_relative_abundance_WT_Rep1,
                                                Intron_relative_abundance_ctr1d_Rep1_per_WT_Rep1 = Intron_relative_abundance_ctr1d_Rep1 / Intron_relative_abundance_WT_Rep1,
                                                Intron_relative_abundance_ctr1d_Rep2_per_WT_Rep1 = Intron_relative_abundance_ctr1d_Rep2 / Intron_relative_abundance_WT_Rep1
                              )
head(ann.median.table)


write.table(ann.median.table, "ann.median.table_Burke.bed", sep = "\t", quote = F, row.names = F)

xy.breaks <- 10^(0:5) 
gg.WT_Rep1.ctrd1.TF43 <- ggplot(ann.median.table, aes(x = WT_Rep1.body, y = ctr1d_Rep1.body)) + geom_point(size = .5) + scale_x_log10(breaks = xy.breaks, expand = c(0.01,0.01), limits = c(1, 10^4.1)) + scale_y_log10(breaks = xy.breaks, expand = c(0.01,0.01), limits = c(1, 10^4.1)) + theme_bw()
gg.WT_Rep1.ctrd1.TF54 <- ggplot(ann.median.table, aes(x = WT_Rep1.body, y = ctr1d_Rep2.body)) + geom_point(size = .5) + scale_x_log10(breaks = xy.breaks, expand = c(0.01,0.01), limits = c(1, 10^4.1)) + scale_y_log10(breaks = xy.breaks, expand = c(0.01,0.01), limits = c(1, 10^4.1)) + theme_bw()
gg.ctrd1.TF43.ctrd1.TF54 <- ggplot(ann.median.table, aes(x = ctr1d_Rep1.body, y = ctr1d_Rep2.body)) + geom_point(size = .5) + scale_x_log10(breaks = xy.breaks, expand = c(0.01,0.01), limits = c(1, 10^4.1)) + scale_y_log10(breaks = xy.breaks, expand = c(0.01,0.01), limits = c(1, 10^4.1)) + theme_bw()


ggsave(gg.WT_Rep1.ctrd1.TF43, filename = "gg.WT_Rep1.ctrd1.TF43_intron_profile.pdf", width = 5, height = 5)
ggsave(gg.WT_Rep1.ctrd1.TF54, filename = "gg.WT_Rep1.ctrd1.TF54_intron_profile.pdf", width = 5, height = 5)
ggsave(gg.ctrd1.TF43.ctrd1.TF54, filename = "gg.ctrd1.TF43.ctrd1.TF54_intron_profile.pdf", width = 5, height = 5)

gg.WT_Rep1.ctr1d.Rep1 <- ggplot(ann.median.table, aes(x = Intron_relative_abundance_WT_Rep1, y = Intron_relative_abundance_ctr1d_Rep1)) + geom_point(size = .1, color = "red") + xlim(0,1.5) + ylim(0, 1.5) + theme_bw()
gg.WT_Rep1.ctr1d.Rep2 <- ggplot(ann.median.table, aes(x = Intron_relative_abundance_WT_Rep1, y = Intron_relative_abundance_ctr1d_Rep2)) + geom_point(size = .1, color = "red") + xlim(0,1.5) + ylim(0, 1.5) + theme_bw()

ggsave(gg.WT_Rep1.ctr1d.Rep1, filename = "gg.WT_Rep1.ctr1d.Rep1_relative_abundance.pdf", width = 5, height = 5)
ggsave(gg.WT_Rep1.ctr1d.Rep2, filename = "gg.WT_Rep1.ctr1d.Rep2_relative_abundance.pdf", width = 5, height = 5)

#########
reliable.indices <- ann.median.table$ctr1d_Rep2.body >= 30 & ann.median.table$ctr1d_Rep1.body >= 30
cols <- ifelse(reliable.indices, "red", "grey")

gg.WT.ctr1d.minus <- ggplot(ann.median.table, aes(x = Intron_relative_abundance_ctr1d_Rep1_minus_WT_Rep1, y = Intron_relative_abundance_ctr1d_Rep2_minus_WT_Rep1)) + geom_point(size = .05, color = cols, alpha = 0.5) + xlim(-1, 1) + ylim(-1, 1) + theme_bw()
gg.WT.ctr1d.per <- ggplot(ann.median.table, aes(x = Intron_relative_abundance_ctr1d_Rep1_per_WT_Rep1,   y = Intron_relative_abundance_ctr1d_Rep2_per_WT_Rep1)) + geom_point(size = .05, color = cols, alpha = 0.5) + xlim(0, 25) + ylim(0, 25) + theme_bw()

ggsave(gg.WT.ctr1d.minus, filename = "gg.WT.ctr1d_relative_abundance_minus.pdf", width = 5, height = 5)
ggsave(gg.WT.ctr1d.per, filename = "gg.WT.ctr1d_relative_abundance_per.pdf", width = 5, height = 5)

cap2_01 <- function(x) {
  capped.x <- ifelse(x>1, 1, ifelse(x < 0, 0, x))
  return(capped.x)
} 

capped.ann.median.table <- transform(ann.median.table, 
                                     Intron_relative_abundance_WT_Rep1 = cap2_01(Intron_relative_abundance_WT_Rep1), 
                                     Intron_relative_abundance_ctr1d_Rep1 = cap2_01(Intron_relative_abundance_ctr1d_Rep1),
                                     Intron_relative_abundance_ctr1d_Rep1_minus_WT_Rep1 = cap2_01(Intron_relative_abundance_ctr1d_Rep1_minus_WT_Rep1), 
                                     Intron_relative_abundance_ctr1d_Rep2_minus_WT_Rep1 = cap2_01(Intron_relative_abundance_ctr1d_Rep2_minus_WT_Rep1)
                                     )


gg.capped.WT.ctr1d.minus <- ggplot(capped.ann.median.table, aes(x = Intron_relative_abundance_ctr1d_Rep1_minus_WT_Rep1, y = Intron_relative_abundance_ctr1d_Rep2_minus_WT_Rep1)) + geom_point(size = .05, color = cols, alpha = 0.5) + xlim(0, 1) + ylim(0, 1) + theme_bw()
ggsave(gg.capped.WT.ctr1d.minus, filename = "gg.capped.WT.ctr1d_relative_abundance_minus.pdf", width = 5, height = 5)

gg.capped.WT_Rep1.ctr1d.Rep1 <- ggplot(capped.ann.median.table, aes(x = Intron_relative_abundance_WT_Rep1, y = Intron_relative_abundance_ctr1d_Rep1)) + geom_point(size = .1, color = "red") + xlim(0,1) + ylim(0, 1) + theme_bw()
gg.capped.WT_Rep1.ctr1d.Rep2 <- ggplot(capped.ann.median.table, aes(x = Intron_relative_abundance_WT_Rep1, y = Intron_relative_abundance_ctr1d_Rep2)) + geom_point(size = .1, color = "red") + xlim(0,1) + ylim(0, 1) + theme_bw()

ggsave(gg.capped.WT_Rep1.ctr1d.Rep1, filename = "gg.capped.WT_Rep1.ctr1d.Rep1_relative_abundance.pdf", width = 5, height = 5)
ggsave(gg.capped.WT_Rep1.ctr1d.Rep2, filename = "gg.capped.WT_Rep1.ctr1d.Rep2_relative_abundance.pdf", width = 5, height = 5)

gg.capped.WT.ctr1d.minus <- ggplot(capped.ann.median.table, aes(x = Intron_relative_abundance_ctr1d_Rep1_minus_WT_Rep1, y = Intron_relative_abundance_ctr1d_Rep2_minus_WT_Rep1)) + geom_point(size = .05, color = cols, alpha = 0.5) + xlim(0, 1) + ylim(0, 1) + theme_bw()
ggsave(gg.capped.WT.ctr1d.minus, filename = "gg.capped.WT.ctr1d_relative_abundance_minus.pdf", width = 5, height = 5)

