###estimate score
install.packages("~/packages/estimate_1.0.13.tar.gz", repos = NULL)
library(estimate)

tpm_gene_log2 <- tpm_gene_log2[!duplicated(tpm_gene_log2$gene),]
rownames(tpm_gene_log2) <- tpm_gene_log2$gene
tpm_gene_log2 <- tpm_gene_log2 %>% select(-c(1,2))
write.table(tpm_gene_log2,file = "estimate_input_exp_mets.txt",col.names = T,row.names = T,sep = "\t",quote = F)

filterCommonGenes(input.f="~/test/immune_score/estimate_input_exp_mets.txt",
                  output.f="~/test/immune_score/estimate_input_exp_mets.gct", id="GeneSymbol")
estimateScore(input.ds = "~/test/immune_score/estimate_input_exp_mets.gct",
              output.ds="~/test/immune_score/estimate_mets_scores.gct",
              platform="illumina")
est <- read.delim("~/test/immune_score/estimate_mets_scores.gct", header=T, skip=2,check.names = F)
est <- as.data.frame(est)
est[,0:5]
rownames(est) <- est$NAME
est$Description <- NULL
est$NAME <- NULL

est=t(est[,3:ncol(est)]) %>% as.data.frame()
est$sample <- rownames(est)
est$sample <- gsub("\\.","-",est$sample)
saveRDS(est,file = "~/test/immune_score/pancancer_estimate_score.rds")