library(dplyr)
library(fgsea)

dt <- readxl::read_xlsx("data/Extrachromosomal DNA is associated with oncogene amplification and poor outcome across multiple cancers.xlsx",sheet = 3)
dt <- dt %>%
  filter(tumor_or_normal!="normal") %>%
  filter(grepl("TCGA",sample_barcode))
dt <- dt %>% group_by(sample_barcode) %>%
  summarise(ecDNA=ifelse(any(sample_classification=="Circular"),"yes","no"))

dt$cancer <- get_cancer_type(dt$sample_barcode)

tcga_trans_counts <- readRDS("~/useful_data/xena_RSEM_expected_counts/tcga_trans_counts.rds")
tcga_trans_counts <- tcga_trans_counts %>%
  select(-sample) %>%
  select(gene,any_of(dt$sample_barcode))

dt <- dt %>% filter(dt$sample_barcode %in% colnames(tcga_trans_counts))
dt %>% group_by(cancer) %>%
  summarise(ecdna_c=sum(ecDNA=="yes"))  %>% filter(ecdna_c > 20) -> summ

need_cancer <- unique(summ$cancer)
library(fgsea)
library(EasyBioinfo)
mhc_i_go <- gmtPathways("data/mhc_i_geneset.gmt")
mhc_ii_go <- gmtPathways("data/mhc_ii_geneset.gmt")
mhc_i_re <- gmtPathways("data/REACTOME_MHC_I.gmt")
mhc_ii_re <- gmtPathways("data/REACTOME_MHC_II.gmt")
pathway <- c(mhc_i_go,mhc_ii_go,mhc_i_re,mhc_ii_re)


res_gsea <- vector("list",7)
names(res_gsea) <- need_cancer

for (i in seq_along(gsea_plot_I_go)){
  
  cancer_dt <- dt %>% filter(cancer == names(gsea_plot_I_go)[i])
  
  gsea_dt <- tcga_trans_counts %>%
    select(gene,any_of(cancer_dt$sample_barcode)) %>%
    as.data.frame()
  gsea_dt <- gsea_dt[!duplicated(gsea_dt$gene),]
  rownames(gsea_dt) <- gsea_dt$gene
  gsea_dt <- gsea_dt[,-1]
  
  cancer_dt <- as.data.frame(cancer_dt)
  ecdna_sample <- cancer_dt %>%  filter(ecDNA=="yes") %>%
    distinct(sample_barcode) %>% filter(sample_barcode %in% colnames(gsea_dt))
  non_ecdna_sample <-  cancer_dt %>% filter(ecDNA=="no") %>%
    distinct(sample_barcode) %>% filter(sample_barcode %in% colnames(gsea_dt))
  gsea_dt <- gsea_dt %>% select(non_ecdna_sample$sample_barcode,ecdna_sample$sample_barcode)
  
  res <- deg_deseq2(gsea_dt,control_label="non_ecDNA",
                    control_counts=nrow(non_ecdna_sample),
                    treatment_lable="ecDNA",treatment_counts=nrow(ecdna_sample),
                    parallel = TRUE,ncores = 18)
  res$gene <- rownames(res)
  ranks <- setNames(res[,"stat"], res[,"gene"])
  fgseaRes <- fgsea::fgsea(pathways = pathway,
                           stats    = ranks,
                           eps      = 1e-10,
                           minSize  = 15,
                           maxSize  = 500)
  res_gsea[[i]] <- fgseaRes
}

saveRDS(res_gsea,file = "data/res_gsea.rds")

