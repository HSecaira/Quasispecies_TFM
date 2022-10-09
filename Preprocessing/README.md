# Pipeline workflow

This directory contains the script `pipeline.sh` for preprocessing raw paired-end reads that will be later used to analyze quasispecies viral populations. 

## Steps 

1. Remove barcodes, primmer sequencues,poor quality bases form the 3'-end of reads, and reads whose length is less than 50 bp using CutAdapt v3.0.
2. Generate quality reports using FastQC v0.11.9 and MultiQC v1.9.
3. Merge forward and reverse reads using Flash2 v2.2.02.
4. Align reads against the bacteriophage $\text{Q}\beta$ reference genome using BWA v0.7.17.
5. Filter reads that contain indels and whose length is less than the region of interest using SAMtools v1.11 and AWK. 
6. Collapse identical reads into a single sequnce keeping a record of their abundance using FASTX-Toolkit v0.0.14. 

![alt text](https://github.com/HSecaira/Quasispecies_TFM/blob/main/Preprocessing/pipelineHorizontal.pdf)
