#!/bin/bash

#SBATCH -c 8

#SBATCH --mem=32G

#SBATCH -t 10:00:00

#SBATCH -J mspc_SOX9N

##### MSPC FOR ID'ING CONSENSUS PEAKS IN ATAC DATASETS. USES .narrowPeak BED6+4 FILES AS INPUT

#Define variables
#OUTPUT_DIR="/proj/jraablab/users/jbrink/sox9/"      # Directory to store mspc output

# Ensure directories exist
#mkdir -p "$OUTPUT_DIR"

#OUTPUT DIR
#OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE_ID}_output"

####run mspc
/proj/jraablab/users/jbrink/sox9/work/mspc/mspc \
  -i "/proj/jraablab/users/jbrink/sox9/aligned_data/peaks/*Sox9P*.narrowPeak" \
  -o "/work/users/j/b/jbrink/processed_atac_sox9p" \
  -r bio \
  -w 1e-6 \
  -s 1e-14 \
  -c 3 \
  -m lowest \
  -d 4 \
