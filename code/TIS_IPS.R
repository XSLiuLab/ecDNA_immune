####################################################
##
##   This R-script can be used to calculate Immunophenoscore (IPS) and generate Immunophenogram from "EXPR.txt" and "IPS_genes.txt"
##   (C) ICBI, Medical University of Innsbruck, Biocenter, Division of Bioinformatics
##   Version 1.0 08.07.2016
##   Needs packages ggplot2,grid,gridExtra
##
####################################################

library(ggplot2)
library(grid)
library(gridExtra)

## calculate Immunophenoscore
ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}
## Read expression data from tab-delimited text file, with official human gene symbols (HGNC) in the first columns
## and expression values (i.e. log2(TPM+1)) for each sample in the other columns
tpm_gene <- readRDS("~/useful_data/xena_RSEM_TPM/tpm_gene.rds")
tpm_gene <- as.data.frame(tpm_gene)
tpm_gene <- tpm_gene[,-1]
tpm_gene <- tpm_gene[!duplicated(tpm_gene$gene),]
rownames(tpm_gene) <- tpm_gene$gene
tpm_gene <- tpm_gene[,-1]

sample_names<-names(tpm_gene)

## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
# For different
IPSG<-read.table("~/test/immune_score/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
unique_ips_genes<-as.vector(unique(IPSG$NAME))

IPS<-NULL
MHC<-NULL
CP<-NULL
EC<-NULL
SC<-NULL
AZ<-NULL

# Gene names in expression file
GVEC<-row.names(tpm_gene)
# Genes names in IPS genes file
VEC<-as.vector(IPSG$GENE)
# Match IPS genes with genes in expression file
ind<-which(is.na(match(VEC,GVEC)))
# List genes missing or differently named
MISSING_GENES<-VEC[ind]
dat<-IPSG[ind,]
if (length(MISSING_GENES)>0) {
  cat("differently named or missing genes: ",MISSING_GENES,"\n")
}
for (x in 1:length(ind)) {
  print(IPSG[ind,])
}

for (i in 1:length(sample_names)) {
  GE<-tpm_gene[[i]]
  mGE<-mean(GE)##所有基因的均值
  sGE<-sd(GE)##所有基因的方差
  Z1<-(tpm_gene[as.vector(IPSG$GENE),i]-mGE)/sGE##IPS gene z score
  W1<-IPSG$WEIGHT
  WEIGHT<-NULL
  MIG<-NULL
  k<-1
  for (gen in unique_ips_genes) {
    MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
    WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
    k<-k+1
  }
  WG<-MIG*WEIGHT
  MHC[i]<-mean(WG[1:10])
  CP[i]<-mean(WG[11:20])
  EC[i]<-mean(WG[21:24])
  SC[i]<-mean(WG[25:26])
  AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
  IPS[i]<-ipsmap(AZ[i])
}
DF<-data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
saveRDS(DF,file = "~/test/immune_score/IPS.rds")

###cal TIS
tpm_gene <- readRDS("~/useful_data/xena_RSEM_TPM/tpm_gene_log2.rds")##GSVA 高斯分布
tpm_gene <- as.data.frame(tpm_gene)
tpm_gene <- tpm_gene[,-2]
tpm_gene <- tpm_gene[!duplicated(tpm_gene$gene),]
rownames(tpm_gene) <- tpm_gene$gene
tpm_gene <- tpm_gene[,-1]

sample_names<-names(tpm_gene)

TIS <- c("CCL5", "CD27", "CD274", "CD276", "CD8A",
         "CMKLR1", "CXCL9", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3",
         "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT","CXCR6")
library(GSVA)
GSVA <- gsva(expr=as.matrix(tpm_gene),
             gset.idx.list=list(TIS), method="gsva",parallel.sz=10)
GSVA <- as.data.frame(t(GSVA))
#save
write.table(GSVA, file="~/test/immune_score/tis_signature.txt", sep="\t", col.names = T, row.names = T, quote = F)
