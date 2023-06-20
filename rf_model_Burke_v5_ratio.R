library(e1071)
library(seqinr)
library(randomForest)
library(party)
library(dplyr)
library(pROC)
library(bedr)

cap2_01 <- function(x) {
  capped.x <- ifelse(x>1, 1, ifelse(x < 0, 0, x))
  return(capped.x)
} 

orig.ann.median.table <- read.table("ann.median.table_Burke.bed", header = T)
merged.ann.median.table <- read.table("orig_Burke_Chen.bed", header = F)

merged.ann.median.table <- transform(merged.ann.median.table, V2 = V54, V3 = V55)[,1:52]
colnames(merged.ann.median.table) <- colnames(orig.ann.median.table)
#full.table <- tt
reliable.indices <- merged.ann.median.table$ctr1d_TF54.body >= 30 & merged.ann.median.table$ctr1d_TF43.body >= 30
tt.Burke <- merged.ann.median.table[reliable.indices,]

#SPAC6C3.03c



#######################
bed.tt <- tt.Burke
bed.tt <- subset(bed.tt, width > 10)
colnames(bed.tt)[1] <- "chr"

levels(bed.tt$chr) <- c("I", "II", "III")
bed.tt[,1] <- as.character(bed.tt[,1])
flanked.bed.tt <- bed.tt
flanked.bed.tt$start <- flanked.bed.tt$start - 9
flanked.bed.tt$end <- flanked.bed.tt$end + 9
#flanked.bed.tt <- flanked.bed.tt[,c(13,6,4,5)]
ff <- get.fasta(bedr.sort.region(flanked.bed.tt[,c(1:6)], check.chr = FALSE), strand = T, fasta = "/data/index/pombase/ASM294v2.fa", check.chr = FALSE, use.name.field = T)

#######################
full.table <- bed.tt
#full.table <- full.table[full.table$WT.upstream > 100 & full.table$WT.downstream > 100,]
full.table <- full.table[order(full.table$Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54),]
full.table$Class <- ifelse(full.table$Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54>0.2, "Sens", "Insens")
cols <- ifelse(full.table$Class == "Sens", "blue", "red")
gg.Sens.minus <- ggplot(full.table, aes(x = Intron_relative_abundance_ctr1d_TF43_minus_WT_TF54,  y = Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54)) + geom_point(size = .05, color = cols, alpha = 0.5) + xlim(-1, 1) + ylim(-1, 1) + theme_bw()
ggsave(gg.Sens.minus, filename = "gg.Sens_minus.pdf", width = 5, height = 5)

capped.full.table <- transform(full.table, Intron_relative_abundance_ctr1d_TF43_minus_WT_TF54 = cap2_01(Intron_relative_abundance_ctr1d_TF43_minus_WT_TF54), 
                                                 Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54 = cap2_01(Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54))
gg.capped.Sens.minus <- ggplot(capped.full.table, aes(x = Intron_relative_abundance_ctr1d_TF43_minus_WT_TF54,  y = Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54)) + geom_point(size = .05, color = cols, alpha = 0.5) + xlim(0, 1) + ylim(0, 1) + theme_bw()
ggsave(gg.capped.Sens.minus, filename = "gg.capped.Sens_minus.pdf", width = 5, height = 5)

#full.table <- rbind(head(full.table, n = 200), tail(full.table, n = 200))
full.table$Class <- as.factor(full.table$Class)

gg.5ss <- ggplot(full.table, aes(X5p.score , Class, fill = Class)) + geom_boxplot() + coord_flip() + theme_bw() #+ xlim(0,10)
ggsave(gg.5ss, filename = "X5p.score_bp_Burke.pdf", width = 7)
gg.3ss <- ggplot(full.table, aes(X3p.score , Class, fill = Class)) + geom_boxplot() + coord_flip() + theme_bw()# + xlim(0,10)
ggsave(gg.3ss, filename = "X3p.score_bp_Burke.pdf", width = 7)

bed.tt <- full.table
levels(bed.tt$chr) <- c("I", "II", "III")
bed.tt[,1] <- as.character(bed.tt[,1])
flanked.bed.tt <- bed.tt
flanked.bed.tt$start <- flanked.bed.tt$start - 9
flanked.bed.tt$end <- flanked.bed.tt$end + 9
#flanked.bed.tt <- flanked.bed.tt[,c(13,6,4,5)]
ff <- get.fasta(bedr.sort.region(flanked.bed.tt[,c(1:6)], check.chr = FALSE), strand = T, fasta = "/data/index/pombase/ASM294v2.fa", check.chr = FALSE, use.name.field = T)


for (class.name in unique(bed.tt$Class)) {
  class.index <- bed.tt$Class == class.name
  print(class.index)
  five.p.res <- count.bases(ff$sequence[class.index], from = 1, to = 20)
  print(five.p.res$count.table)
  five.p.bits.logo <- ggseqlogo(five.p.res$count.table,  method = 'bits'); ggsave(five.p.bits.logo, height = 3, width = 12, filename = paste0(class.name, "_Burke_ann_5pSS_logo_bits.pdf"))
  three.p.res <- count.bases(ff$sequence[class.index], from = -21, to = NULL)
  three.p.bits.logo <- ggseqlogo(three.p.res$count.table,  method = 'bits'); ggsave(three.p.bits.logo, height = 3, width = 12, filename = paste0(class.name,"_Burke_ann_3pSS_logo_bits.pdf"))
  
  
  cc <- cbind(bed.tt[class.index,], seq = ff$sequence[class.index])
  cc <- cc[order(cc$Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54, decreasing = T),]
  write.table(cc, paste0(class.name, "_n", nrow(cc),"_bed_tt.tsv"), sep = "\t", quote = F, row.names = F)
}


######
Insens.tt <- subset(bed.tt, Class == "Insens")
ordered.Insens.tt <- Insens.tt[order(Insens.tt$WT_TF54.upstream, decreasing = T),]
print(nrow(Insens.tt))
  bed.tt <- ordered.Insens.tt
  levels(bed.tt$chr) <- c("I", "II", "III")
  bed.tt[,1] <- as.character(bed.tt[,1])
  flanked.bed.tt <- bed.tt
  flanked.bed.tt$start <- flanked.bed.tt$start - 9
  flanked.bed.tt$end <- flanked.bed.tt$end + 9
  upscaled.bed.tt <- NULL
    for (i in 1:nrow(flanked.bed.tt)) {
      #full.table$Intron_cut_abundance_ctr1d_minus_WT
      print(i)
      if (flanked.bed.tt$Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54[i] <= 0) next
      scale.value <- round(100*flanked.bed.tt$Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54[i]*flanked.bed.tt$WT_TF54.upstream[i])
      print(scale.value)
      # upscaled.bed.tt <- rbind(upscaled.bed.tt, flanked.bed.tt[i,c(1:6)][rep(1,round(full.table$Intron_cut_abundance_ctr1d_minus_WT[i]*flanked.bed.tt$WT.upstream[i])),])
      #upscaled.bed.tt <- rbind(upscaled.bed.tt, flanked.bed.tt[i,c(1:6)][rep(1,round(flanked.bed.tt$WT.upstream[i])),])
      upscaled.bed.tt <- rbind(upscaled.bed.tt, flanked.bed.tt[i,c(1:6)][rep(1,scale.value),])
    }
    up.scaled.ff <- get.fasta(bedr.sort.region(upscaled.bed.tt, check.chr = FALSE), strand = T, fasta = "/data/index/pombase/ASM294v2.fa", check.chr = FALSE, use.name.field = T)

  aa <- summary(as.factor(unlist(lapply(up.scaled.ff$sequence, substring, 1,11))), maxsum = 100000)
  aaa <- sort(aa, decreasing = T)
  ww <- wordcloud(names(aaa)[1:100], aaa[1:100], scale = c(5, 0.25))
  
    
  five.p.res <- count.bases(up.scaled.ff$sequence, from = 1, to = 20)
  print(five.p.res$count.table)
  five.p.bits.logo <- ggseqlogo(five.p.res$count.table,  method = 'bits'); ggsave(five.p.bits.logo, height = 3, width = 12, filename = paste0("Insense_scaled_to_Exp_Burke_ann_5pSS_logo_bits.pdf"))
  three.p.res <- count.bases(up.scaled.ff$sequence, from = -21, to = NULL)
  three.p.bits.logo <- ggseqlogo(three.p.res$count.table,  method = 'bits'); ggsave(three.p.bits.logo, height = 3, width = 12, filename = paste0("Insense_scaled_to_Exp_Burke_ann_3pSS_logo_bits.pdf"))
  
  sorted.flanked.bed.tt <- bedr.sort.region(flanked.bed.tt[,1:6], check.chr = FALSE)
  sorted.ff <- get.fasta(sorted.flanked.bed.tt, strand = T, fasta = "/data/index/pombase/ASM294v2.fa", check.chr = FALSE, use.name.field = T)
  
  mm <- cbind(bed.tt, seq = sorted.ff$sequence)
  mm <- mm[order(mm$Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54, decreasing = T),]
  sorted.Insens.tt <- cbind(flanked.bed.tt[rownames(sorted.flanked.bed.tt),], seq = sorted.ff$sequence)
  write.table(sorted.Insens.tt, paste0("Insense_", nrow(mm), "_bed_tt.tsv"), sep = "\t", quote = F, row.names = F)
  write.table(Insens.tt, paste0("Insense_", nrow(mm), "_orig.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
  
  
#######
  ######
  bed.tt <- full.table
  Sens.tt <- subset(bed.tt, Class == "Sens")
  ordered.Sens.tt <- Sens.tt[order(Sens.tt$WT_TF54.upstream, decreasing = T),]
  print(nrow(Sens.tt))
  bed.tt <- ordered.Sens.tt
  levels(bed.tt$chr) <- c("I", "II", "III")
  bed.tt[,1] <- as.character(bed.tt[,1])
  flanked.bed.tt <- bed.tt
  flanked.bed.tt$start <- flanked.bed.tt$start - 9
  flanked.bed.tt$end <- flanked.bed.tt$end + 9
  upscaled.bed.tt <- NULL
  for (i in 1:nrow(flanked.bed.tt)) {
    #full.table$Intron_cut_abundance_ctr1d_minus_WT
    print(i)
    if (flanked.bed.tt$Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54[i] <= 0) next
    scale.value <- round(100*flanked.bed.tt$Intron_relative_abundance_ctr1d_TF54_minus_WT_TF54[i]*flanked.bed.tt$WT_TF54.upstream[i])
    print(scale.value)
    # upscaled.bed.tt <- rbind(upscaled.bed.tt, flanked.bed.tt[i,c(1:6)][rep(1,round(full.table$Intron_cut_abundance_ctr1d_minus_WT[i]*flanked.bed.tt$WT.upstream[i])),])
    #upscaled.bed.tt <- rbind(upscaled.bed.tt, flanked.bed.tt[i,c(1:6)][rep(1,round(flanked.bed.tt$WT.upstream[i])),])
    upscaled.bed.tt <- rbind(upscaled.bed.tt, flanked.bed.tt[i,c(1:6)][rep(1,scale.value),])
  }
  up.scaled.ff <- get.fasta(bedr.sort.region(upscaled.bed.tt, check.chr = FALSE), strand = T, fasta = "/data/index/pombase/ASM294v2.fa", check.chr = FALSE, use.name.field = T)
  
  five.p.res <- count.bases(up.scaled.ff$sequence, from = 1, to = 20)
  print(five.p.res$count.table)
  five.p.bits.logo <- ggseqlogo(five.p.res$count.table,  method = 'bits'); ggsave(five.p.bits.logo, height = 3, width = 12, filename = paste0("Sense_scaled_to_Exp_Burke_ann_5pSS_logo_bits.pdf"))
  three.p.res <- count.bases(up.scaled.ff$sequence, from = -21, to = NULL)
  three.p.bits.logo <- ggseqlogo(three.p.res$count.table,  method = 'bits'); ggsave(three.p.bits.logo, height = 3, width = 12, filename = paste0("Sense_scaled_to_Exp_Burke_ann_3pSS_logo_bits.pdf"))

  sorted.flanked.bed.tt <- bedr.sort.region(flanked.bed.tt[,1:6], check.chr = FALSE)
  sorted.ff <- get.fasta(sorted.flanked.bed.tt, strand = T, fasta = "/data/index/pombase/ASM294v2.fa", check.chr = FALSE, use.name.field = T)

  ###
  library(wordcloud)
  
  options(max.print = 1000000)
  aa <- summary(as.factor(unlist(lapply(up.scaled.ff$sequence, substring, 1,11))), maxsum = 100000)
  aaa <- sort(aa, decreasing = T)
  ww <- wordcloud(names(aaa)[1:100], aaa[1:100], scale = c(5, 0.25))
  ##
  
  #ordered.Sens.tt$index <- paste0(ordered.Sens.tt$gene_id, "::", ordered.Sens.tt$chr, ":", ordered.Sens.tt$start, "-", ordered.Sens.tt$end, "(", ordered.Sens.tt$strand, ")")
  #merged.Sens.tt <- merge(ordered.Sens.tt, as.data.frame(sorted.ff), by = "index")
  sorted.Sens.tt <- cbind(flanked.bed.tt[rownames(sorted.flanked.bed.tt),], seq = sorted.ff$sequence)
  write.table(sorted.Sens.tt, paste0("Sense_", nrow(cc), "_bed_tt.tsv"), sep = "\t", quote = F, row.names = F)
  
  write.table(Sens.tt, paste0("Sense_", nrow(cc), "_orig.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
  
  
