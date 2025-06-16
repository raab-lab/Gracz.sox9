library(stringr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

#load vsd, granges and results RDS

anno <- annotatePeak(granges,
                     TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                     annoDb = "org.Mm.eg.db",
                     tssRegion = c(-3000, 3000)
)
meta <- mcols(anno@anno)
cts_mat <- assay(vsd)
cts_meta <- cbind(meta, cts_mat)
rownames(cts_meta) <- 1:nrow(cts_meta)
rownames(res_shk) <- 1:nrow(res_shk)
cts_meta <- cts_meta[!grepl("^Gm|^Rik|^Mup", cts_meta$SYMBOL), ]

to_heat <- cts_meta[order(res_filt$padj)[1:3000], ]
to_heat <- cts_meta[order(res_filt$log2FoldChange), ]

res_filt <- res_shk[rownames(cts_meta), ]


res_filt$SYMBOL <- cts_meta$SYMBOL
saveRDS(cts_meta,"/Users/jbrink/Gracz.sox9/src/atac/inj_atac_ctsmeta.rds")
saveRDS(res_filt, "/Users/jbrink/Gracz.sox9/src/atac/inj_atac_resFilt.rds")