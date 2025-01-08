# FFT-trial-DK
This is the repository for R scripts used in the analysis of 16S rRNA and shotgun microbiome data analysis generated in the Danish FFT trial

*updated 08/01/2025, Mattia Pirolo*

## Objective
To examine the impact of fecal filtrate transplantation (FFT) on prevention of post-weaning diarreha in piglets.

## Study design
![FFT trial x paper](https://github.com/user-attachments/assets/5e6c3287-403a-40c6-8bb7-f340bb63b65b)

## Repository organization
### Data
This folder contains all datafiles:
- **ps_16S.rds**: phyloseq object in R data format for the 16S rRNA data analysis
- **ps_MAG.rds**: phyloseq object in R data format for the MAG data analysis
### R script
This folder contains the R script:
- **FFT_study_analysis.R**: R script for visualization of taxonomic composition, alpha- and beta-diversity analysis, and differential abundance analysis using [LEfSe](https://www.bioconductor.org/packages/release/bioc/html/lefser.html)
