library(plyranges)
library(tidyverse)
library(DESeq2)
library(csaw)
library(stringr)

union_peaks <- read_bed(file.path('/work/users/j/b/jbrink/union_uninj_inj_atac.bed'))

bams <- list.files(
  path = '/proj/jraablab/users/jbrink/sox9/aligned_data/bams/filtered',
  full.names = TRUE
)

bam_files <- bams[grep("Ms1.*(chow|ddc).*\\.bam$", bams)]

# CSAW Setup
param = readParam(pe = 'both', 
                  dedup = F, 
                  minq = 10) # set this to 10 which matches my bam parameters, higher and we lose some numbers

peak_counts <- regionCounts(bams_exp, 
                            union_peaks, 
                            param = param
)

peak_counts

#checking dims, etc + counts distro b/c not sure what to expect with peak counts
assayNames(peak_counts)

dim(assay(peak_counts))

hist(assay(peak_counts), 
     breaks=5, 
     main="Distribution of Counts", 
     xlab="Counts"
)

count_matrix <- assay(peak_counts)

# Create a data frame with filenames and their base names
ss <- data.frame(filename = bams_exp, 
                 names = basename(bams_exp)
)
print(ss)

# Mutate and separate the short_name column
ss <- ss %>%
  mutate(short_name = str_replace(names, '.aligned.bam', '')) %>%
  separate(short_name, into = c('airtable_id',
                                'Sample', 
                                'cell_line', 
                                'condition', 
                                'rep'), 
                                sep = '_'
  )

# Print the resulting df
print(ss)

# Join the `ss` data frame with the column data of `peak_norm`
colData(peak_counts) <- DataFrame(left_join(as.data.frame(colData(peak_counts)), 
                                          ss, 
                                          by = c('bam.files' = 'filename'))
)

# Normalize counts using DESeq2 
dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                              colData = colData(peak_counts), 
                              design = ~condition)
dds$condition <- relevel(dds$condition, ref = 'chow')

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

des <-DESeq(dds)

res <- results(des)
vsd <- vst(dds)

plotPCA(vsd, 
        intgroup = "condition")

save(des, file = '/work/users/j/b/jbrink/atac_uninj_inj_de.Rda')
