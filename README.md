# FFT-trial-DK
This is the repository for R scripts used in the analysis of 16S rRNA and shotgun microbiome data analysis generated in the Danish FFT trial

*updated 08/01/2025, Mattia Pirolo*

## Objective
To examine the impact of fecal filtrate transplantation (FFT) on prevention of post-weaning diarreha in piglets.

## Study design
![image](https://github.com/mpirolo/Pilot-repository-16S/assets/54710620/c3fb1634-d19b-4120-b716-2d2f963cc0e7)

## Repository organization
### Data
This folder contains all datafiles:
- **ps_FMT.rds**: phyloseq object in R data format for the analysis
- **metadata.tsv**: tab-separated values file with relevant metadata for the analysis
### R scripts
This folder contains all R scripts:
- **composition.R**: R script for visualization of taxonomic composition
- **diversity.R**: R script for alpha- and beta-diversity analysis
- **DESeq2.R**: R script for differential abundance analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
