```{r}
library(tidyverse)
library(plyranges)
library(Biostrings)
```

```{r union peak calls}
library(rtracklayer)

#read in bed file as dataframe of consensus peaks from each conditon
uninj_consensus <- read.csv("/Users/jbrink/Desktop/Lab2024/gracz.sox9_archive/atac/uninj_ConsensusPeaks.bed", header = TRUE, sep="\t")
inj_consensus <- read.csv("/Users/jbrink/Desktop/Lab2024/gracz.sox9_archive/atac/inj_ConsensusPeaks.bed", header = TRUE, sep="\t")


#convert df to bed for plyranges 
uninj_consen_peak <- makeGRangesFromDataFrame(uninj_consensus,
                         keep.extra.columns=TRUE,
                         ignore.strand = FALSE,
                         na.rm = TRUE)

inj_consen_peak <- makeGRangesFromDataFrame(inj_consensus,
                         keep.extra.columns=TRUE,
                         ignore.strand = FALSE,
                         na.rm = TRUE)

# make a union of the consensus peaks
all_peaks <- bind_ranges(uninj_consen_peak, inj_consen_peak) |> reduce_ranges()

# write out union peak set.
write_bed(all_peaks, file.path("/Users/jbrink/Desktop/Lab2024/gracz.sox9_archive/atac", 'union_uninj_inj_atac.bed'))
```

