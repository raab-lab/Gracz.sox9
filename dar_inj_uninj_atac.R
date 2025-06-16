library(tidyverse)
library(DESeq2)
library(stringr)

#This script performs differential abundance analysis on peaks called in ATACseq injured and uninjured datasets. 

load("~/Downloads/granges.RData")
load("~/Downloads/peak_counts.RData")
se <- peak_counts

deseq_from_se <- function(se, condition_vec, batch_vec){
  
  count_mat <- assay(se)
  
  bam_files <- colData(se)$bam.files
  sample_ids <- gsub(".*_(Ms[0-9]+).*", "\\1", bam_files)
  
  colnames(count_mat) <- sample_ids
  
  sample_info <- data.frame(
    row.names = sample_ids,
    condition = condition_vec,
    batch = batch_vec
  )
  
  #stop if rownames do not align with count matrix
  stopifnot(all(rownames(sample_info) == colnames(count_mat)))
  
  #make deseqdataset
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = sample_info,
    design = ~batch + condition
  )
  
  return(list(count_mat = count_mat,
              sample_info = sample_info,
              dds = dds
  ))
}

conditions <- c("injured","injured","injured", "uninjured", "injured", "uninjured", "uninjured", "uninjured")
batches <- c("B","B", "A", "A", "A", "B", "A", "A")

result <- deseq_from_se(se, conditions, batches)
dds <- result$dds

dds$condition <- relevel(dds_atac$condition, ref = "uninjured")
des <- DESeq(dds)
vsd <- vst(des)
res <- results(des)
summary(res)

res_shk <- lfcShrink(des, 
                     coef = "condition_injured_vs_uninjured", 
                     type = "apeglm"
)

#save objects as RDS and csvs
saveRDS(res, "/Users/jbrink/Gracz.sox9/src/atac/inj_atac_res.rds")
saveRDS(vsd, "/Users/jbrink/Gracz.sox9/src/atac/inj_atac_vsd.rds")
saveRDS(res_shk, "/Users/jbrink/Gracz.sox9/src/atac/inj_atac_resShk.rds")

write.csv(res, "/Users/jbrink/Gracz.sox9/src/atac/inj_atac_res.csv")
write.csv(res_shk, "/Users/jbrink/Gracz.sox9/src/atac/inj_atac_resShk.csv")

