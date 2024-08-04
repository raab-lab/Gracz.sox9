# RNA seq for damaged and undamaged (DDC vs Chow) mouse sox9 cells

suppressPackageStartupMessages({
  library(tidyverse) 
  library(DESeq2) # main package for differential expression
  # Helper functions for reading in count data 
  library(limma)
  library(tximeta)
  library(janitor)
})


# Lazy but do this for the time being
input_dir <- 'data/source_data/rna/quants/'
output_dir <- 'data/processed_data/rna/' # for saving the DES objects
final_dir <- 'data/derived_data/rna/'

for (i in c(output_dir, final_dir)) { 
  if(!dir.exists(i))  {dir.create(i, recursive = T)}
}


# Import 
cda <- read_csv('rna_samples.csv') |> clean_names()
cda$names <- paste('1',cda$sample_id, cda$cell_line, cda$treatment, cda$replicate, sep = '_' )
cda$files <- file.path(input_dir, cda$names, '/quant.sf'  ) 

# use tximeta to import salmon data
txi <- tximeta(coldata = cda, type = 'salmon')

# sumamrise Tx to Genes
se <- summarizeToGene(txi)
dds <- DESeqDataSet(se, design =  ~   treatment)
des <- DESeq(dds)
res <- lfcShrink(des, coef = 2, type = 'apeglm', format = 'GRanges', saveCols=2)
table(res$padj < 0.05, res$log2FoldChange> 0)

###############################################################################
# QC plots
vst <- varianceStabilizingTransformation(dds, blind = T)
plotPCA(vst, intgroup = c('treatment', 'replicate', 'sample_id'), returnData = T) |> 
  ggplot(aes(x = PC1, y = PC2, color = treatment, shape = as.factor(replicate))) + geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = sample_id))

################################################################################
# Save the rda from the deseq datasets
save(des_hlf_ash, des_hlf_men, des_plc, file = 'data/processed_data/rna/des_obj_jr.Rda')

write_tsv(plc_ko_v_wt |> as.data.frame() |> rownames_to_column(), file = 'data/processed_data/rna/plc_ko_v_wt.tsv')

################################################################################
rowRanges(res)

