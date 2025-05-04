library(DiffBind)


# bam file paths
AR_bams = list.files("/Users/prakashsah/stelloo_fig2/data/AR", pattern="^SRR.*\\.sorted.bam$", full.names=TRUE)
H3K27ac_bams = list.files("/Users/prakashsah/stelloo_fig2/data/H3K27ac", pattern="^SRR.*\\.sorted.bam$", full.names=TRUE)
H3K4me3_bams = list.files("/Users/prakashsah/stelloo_fig2/data/H3K4me3", pattern="^SRR.*\\.sorted.bam$", full.names=TRUE)
H3K27me3_bams = list.files("/Users/prakashsah/stelloo_fig2/data/H3K27me3", pattern="^SRR.*\\.sorted.bam$", full.names=TRUE)

# peak file paths
AR_peaks = list.files("/Users/prakashsah/stelloo_fig2/data/AR/peak", pattern = "*.narrowPeak$", full.names = TRUE)
H3K27ac_peaks = list.files("/Users/prakashsah/stelloo_fig2/data/H3K27ac/peak", pattern = "*.narrowPeak$", full.names = TRUE)
H3K4me3_peaks = list.files("/Users/prakashsah/stelloo_fig2/data/H3K4me3/broadpeak", pattern = "*.broadPeak$", full.names = TRUE)
H3K27me3_peaks = list.files("/Users/prakashsah/stelloo_fig2/data/H3K27me3/peak", pattern = "*.narrowPeak$", full.names = TRUE)

# create sample sheet data frame for diffbind
samples <- data.frame(
  SampleID = c(
    sub("\\.sorted\\.bam$", "", basename(AR_bams)),
    sub("\\.sorted\\.bam$", "", basename(H3K27ac_bams)),
    sub("\\.sorted\\.bam$", "", basename(H3K4me3_bams)),
    sub("\\.sorted\\.bam$", "", basename(H3K27me3_bams))
  ),
  Condition = rep(c("Case", "Case", "Control", "Control"), 4),
  Factor = rep(c("AR", "H3K27ac", "H3K4me3", "H3K27me3"), each = 4),
  bamReads = c(AR_bams, H3K27ac_bams, H3K4me3_bams, H3K27me3_bams),
  Peaks = c(AR_peaks, H3K27ac_peaks, H3K4me3_peaks, H3K27me3_peaks),
  PeakCaller = rep("macs", 16)
)

# Read in data and create diffbind object
diff_obj = dba(sampleSheet = samples)
dba.plotHeatmap(diff_obj, attributes = DBA_FACTOR) # heatmap based on peak occupancy
dba.plotPCA(diff_obj, DBA_FACTOR) # PCA plot based on peak occupancy

# create count object and plot heatmap and PCA
diff_obj_count = dba.count(diff_obj)
dba.plotHeatmap(diff_obj_count, attributes = DBA_FACTOR)
dba.plotPCA(diff_obj_count, DBA_FACTOR)




