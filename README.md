# ChIP-seq Analysis Challenge: Subset of Stelloo et al. (2018)

This project is a simplified ChIP-seq analysis using a subset of Stelloo et al. 2018 (GSE120738) to run locally on a Mac. This was done as part of the challenge from Ming Tommy Tang’s CHIPSeq data analysis tutorial.

---

## 🔍 Overview

- **Dataset**: Stelloo et al. 2018 ([GSE120738](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120738))
- **ChIP Factors**: AR, H3K4me3, H3K27ac, H3K27me3
- **Samples**: 2 cases, 2 controls, 2 inputs per factor (patient-matched)
- **Raw data**: Downloaded fastq files fom from ENA (metadata from SRArun selector)

---

## 🛠️ Workflow

### Preprocessing (Bash)
- Fastq download, QC (`FastQC`), alignment (`BWA`), bigWig generation
- Peak calling using `MACS3` (broad peaks for H3K4me3; peaks called against merged input for each factor)

### Analysis (RStudio)
- **Fig 2a**: IGV used to visualize bigWig files at genomic regions specified in the paper, saved as images
- **Fig 2b–c**: Peak lengths extracted from MACS3 output using `awk`; merged and visualized via density plots in R. AR peaks are narrower as expected; differences in 2c likely due to fewer samples and peak caller.
- **Fig 2d–e**: `DiffBind` used for PCA and heatmap; input = sample sheet with BAM and peak file paths
- **Fig 2g**: `ChIPseeker` used to annotate consensus peaksets generated using DiffBind
- **Fig 2h**: Top motifs identified using Galaxy SeqPos (JASPAR database), then visualized in R

---

## 📦 Tools Used

- Bash, BWA, FastQC, MACS3
- R: DiffBind, ChIPseeker, ggplot2
- IGV (Web/Desktop)
- Galaxy SeqPos

---


