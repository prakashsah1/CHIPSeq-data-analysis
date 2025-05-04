library(ChIPseeker)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)

#create consensus peak file for each CHIP factor
##subset diff_obj object by CHIP factor
AR_obj <- dba(diff_obj, mask = diff_obj$samples$Factor == "AR")
H3K27ac_obj <- dba(diff_obj, mask = diff_obj$samples$Factor == "H3K27ac")
H3K4me3_obj <- dba(diff_obj, mask = diff_obj$samples$Factor == "H3K4me3")
H3K27me3_obj <- dba(diff_obj, mask = diff_obj$samples$Factor == "H3K27me3")

##generate consensus peakset for each CHIP factor (add to the AR_obj)
AR_obj <- dba.peakset(AR_obj, consensus = DBA_FACTOR, minOverlap = 2)
H3K27ac_obj <- dba.peakset(H3K27ac_obj, consensus = DBA_FACTOR, minOverlap = 2)
H3K4me3_obj <- dba.peakset(H3K4me3_obj, consensus = DBA_FACTOR, minOverlap = 2)
H3K27me3_obj <- dba.peakset(H3K27me3_obj, consensus = DBA_FACTOR, minOverlap = 2)

##retrieve consensus peak files for each CHIP factor
AR_consensus <- dba.peakset(AR_obj, bRetrieve = TRUE)
H3K27ac_consensus <- dba.peakset(H3K27ac_obj, bRetrieve = TRUE)
H3K4me3_consensus <- dba.peakset(H3K4me3_obj, bRetrieve = TRUE)
H3K27me3_consensus <- dba.peakset(H3K27me3_obj, bRetrieve = TRUE)

#annotate peaks 
AR_annoPeak = annotatePeak(AR_consensus, tssRegion = c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db")
H3K27ac_annoPeak = annotatePeak(H3K27ac_consensus, tssRegion = c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db")
H3K4me3_annoPeak = annotatePeak(H3K4me3_consensus, tssRegion = c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db")
H3K27me3_annoPeak = annotatePeak(H3K27me3_consensus, tssRegion = c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db")

# plot piechart for annotated peaks
plotAnnoPie(AR_annoPeak)
plotAnnoPie(H3K27ac_annoPeak)
plotAnnoPie(H3K4me3_annoPeak)
plotAnnoPie(H3K27me3_annoPeak)
