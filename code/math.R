setwd("e:/ecDNA/MAF/")
library(mclust)
library(sigminer)
library(NMF)
library(stringr)
library(maftools)
a<-list.dirs()
a<-a[-1]
MATH_Result <- data.frame(NA,NA)
colnames(MATH_Result) <- c("d","sample_id")

for(i in 1:length(a)){ 
  a2<-list.files(a[i])
  b<-grep("maf.gz$",a2)
  a2<-a2[b]
  laml.maf = str_c(a[i],a2,sep="/")
  laml <- read.maf(maf = laml.maf,isTCGA=T)
  laml@data$VAF <- laml@data$t_alt_count/laml@data$t_depth 
  d <- matrix()
  sample_id <- names(table(laml@data$Tumor_Sample_Barcode))
  sample_id <- setdiff(sample_id,names(which(table(laml@data$Tumor_Sample_Barcode) < 4)))
  for (sample in sample_id) {
    res <- inferHeterogeneity(maf = laml, tsb = sample, vafCol = 'VAF') 
    d[sample] <- cbind(unique(res$clusterData$MATH),names(table(res$clusterData$Tumor_Sample_Barcode)))
  }
  d <- as.data.frame(d)
  d$sample_id <- rownames(d)
  MATH_Result <- rbind(MATH_Result,d)
}
MATH <- MATH_Result
MATH <- na.omit(MATH)
saveRDS(MATH,file = "e:/project/ecDNA/MAF_Figure/TCGA_MATH_Value.rds")