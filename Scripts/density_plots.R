library(ggplot2)
library(dplyr)

# Read combined peaks sizes for all CHIP factors
AR <- read.table("/Users/prakashsah/stelloo_fig2/data/AR/peak/AR_all_peaks.txt", col.names = "Size")
H3K27ac <- read.table("/Users/prakashsah/stelloo_fig2/data/H3K27ac/peak/H3K27ac_all_peaks.txt", col.names = "Size")
H3K4me3 <- read.table("/Users/prakashsah/stelloo_fig2/data/H3K4me3/broadpeak/H3K4me3_all_broadpeaks.txt", col.names = "Size")
H3K27me3 <- read.table("/Users/prakashsah/stelloo_fig2/data/H3K27ac/peak/H3K27ac_all_peaks.txt", col.names = "Size")

# add peak identifier
AR$Group <- "AR"
H3K27ac$Group <- "H3K27ac"
H3K4me3$Group <- "H3K4me3"
H3K27me3$Group <- "H3K27me3"

# combine AR, H3K27ac and H3K4me3 peaks data
combined_peaks <- bind_rows(AR, H3K27ac, H3K4me3)

# Plot combined peaks
ggplot(combined_peaks, aes(x = Size, color = Group)) +
  geom_density(alpha = 0.4) +
  scale_color_manual(values = c("darkgreen", "darkorange", "deeppink")) +
  labs(title = "Peak Size Density Plot", x = "Peak Size (bp)", y = "Density") +
  coord_cartesian(xlim = c(0, 3000)) +
  theme_minimal()

# Plot H3K27me3 peak size
ggplot(H3K27me3, aes(x = Size)) +
  geom_density(color = "steelblue", size = 1) + 
  labs(x = "Peak Size (bp)", y = "Density") +
  theme_minimal()
