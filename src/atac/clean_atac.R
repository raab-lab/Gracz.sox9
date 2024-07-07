```{r}
# Make consensus atac peaks using mspc
# must have dotnet available
# from console
# module load dotnet 
#make sure you are using bioconductor 3.18 or higher and rsmpc 1.8
library(rmspc)
library(tidyverse)
library(plyranges)
library(Biostrings)

#Sys.setenv(PATH = paste("/usr/local/share/dotnet/dotnet", Sys.getenv("PATH"), sep=""))
```

```{r}
# Output location
output_dir <- 'Users/jbrink/Desktop/Lab2024/gracz.sox9_archive/peak_sets/atac_consensus/'
if(!dir.exists(output_dir) ) { dir.create(output_dir, recursive = T)} 

# Include an excluded region list
col_names <- c('seqnames', 'start', 'end', 'name', 'score', 'strand') 

# I'm excluding these anyway
exclude_regions <- read_tsv('/Users/jbrink/Desktop/Lab2024/gracz.sox9_archive/atac/mm10-exclusion.v2.bed', col_names = col_names) |> as_granges()


# Takes pattern for identifying a group of peaks 
# Returns a named list of peaks for that antibody
process_peak_group <- function(pattern) {
  files <- list.files(path = 'Users/jbrink/Desktop/Lab2024/gracz.sox9_archive/atac/peaks/', pattern = pattern, full.names = T)  
  names <- basename(files)
  names <- lapply(names, function(x)  (unlist(str_split(x, pattern = '_') )  ) ) 
  names <- sapply(names, function(x) paste(x[3:6], collapse = '_') )
  peaks <- lapply(files, function(x) read_narrowpeaks(x) |> filter_by_non_overlaps(exclude_regions, maxgap=500) )  
  names(peaks) <- names
  return(peaks) 
}
```


```{r}
# Get list of peak calls for each condition
uninj_peaks <- process_peak_group('.*chow.*.narrowPeak')
ddc_peaks <- process_peak_group('.*DDC.*.narrowPeak')

# Create consensus set for each antibody_sample ##############################
# need at least 2 replicates for this to make sense 

create_consensus <- function(peak_list, c) {
  peak_consensus <- mspc(input = GRangesList(peak_list), 
                         replicateType = 'Bio', 
                         stringencyThreshold = 1e-14,
                         weakThreshold = 1e-6, 
                         keep = FALSE, 
                         GRanges = TRUE,
                         multipleIntersections = "Lowest",
                         c = c, 
                         alpha = 0.05
  )
  return(peak_consensus$GRangesObjects$ConsensusPeaks)
}

uninj_consensus <- create_consensus(uninj_peaks, c = 4)
ddc_consensus <- create_consensus(ddc_peaks, c = 4)

# make a union of the consensus peaks
all_peaks <- bind_ranges(dmso_consensus, sndx_consensus) |> reduce_ranges()

# write out peak sets
write_bed(uninj_consensus, file.path(output_dir, 'uninj_consensus_atac.bed'))
write_bed(ddc_consensus, file.path(output_dir, 'ddc_consensus_atac.bed'))
write_bed(all_peaks, file.path(output_dir, 'union_atac.bed'))
```

