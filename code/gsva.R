library(dplyr)
library(fgsea)
library(GSVA)

tpm_gene <- readRDS("~/useful_data/xena_RSEM_TPM/tpm_gene_log2.rds")
tpm_gene <- as.data.frame(tpm_gene)
tpm_gene <- tpm_gene %>% select(-id)
tpm_gene <- tpm_gene[!duplicated(tpm_gene$gene),]
rownames(tpm_gene) <- tpm_gene$gene
tpm_gene <- tpm_gene[,-1]

mhc_i_go <- gmtPathways("data/mhc_i_geneset.gmt")
mhc_ii_go <- gmtPathways("data/mhc_ii_geneset.gmt")
mhc_i_re <- gmtPathways("data/REACTOME_MHC_I.gmt")
mhc_ii_re <- gmtPathways("data/REACTOME_MHC_II.gmt")
pathway <- c(mhc_i_go,mhc_ii_go,mhc_i_re,mhc_ii_re)

sapply(1:length(pathway), function(i)
  paste0(names(pathway)[[i]], " ",
         sum(pathway[[i]] %in% rownames(tpm_gene)), "/",
         length(pathway[[i]])))

dt <- readxl::read_xlsx("data/Extrachromosomal DNA is associated with oncogene amplification and poor outcome across multiple cancers.xlsx",sheet = 3)
dt <- dt %>%
  filter(tumor_or_normal!="normal") %>%
  filter(grepl("TCGA",sample_barcode))
dt <- dt %>% group_by(sample_barcode) %>%
  summarise(ecDNA=ifelse(any(sample_classification=="Circular"),"yes","no"))
library(NeoEnrichment)
dt$cancer <- get_cancer_type(dt$sample_barcode)

library(parallel)
tpm_gene <- tpm_gene %>%
  select(any_of(dt$sample_barcode))
dt <- dt %>% filter(sample_barcode %in% colnames(tpm_gene))
dt %>% group_by(cancer) %>% summarise(c=sum(ecDNA=="yes")) %>%
  filter(c>20)-> summ
dt <- dt %>% filter(cancer %in% summ$cancer)
tpm_gene <- tpm_gene %>%
  select(any_of(dt$sample_barcode))

GSVA <- gsva(expr=as.matrix(tpm_gene), gset.idx.list=pathway,
             verbose=FALSE, mx.diff=1, method="gsva",parallel.sz=15)
saveRDS(GSVA,file = "data/GSVA_pathway.rds")

