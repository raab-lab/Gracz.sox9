library(plyranges)
library(tidyverse)
library(DESeq2)
library(csaw)
library(stringr)

union_peaks <- read_bed(file.path('/work/users/j/b/jbrink/union_uninj_inj_atac.bed'))

bams <- list.files(
  path = '/proj/jraablab/users/jbrink/sox9/aligned_data/bams/filtered/',
  full.names = TRUE
)
bams_exp <- bams[grep("Ms1.*(chow|ddc).*\\.bam$", bams)]

# CSAW Setup
param = readParam(pe = 'both', dedup = F, minq = 10) # set this to 10 which matches my bam parameters, higher and we lose some numbers
peak_counts <- regionCounts(bams_exp, union_peaks, param = param)
