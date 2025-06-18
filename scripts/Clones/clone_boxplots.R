# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# ---------------------------
# 1. Set up sample & working directory
# ---------------------------
sample_id <- "DU34"  # ← change for each sample
setwd("~/Desktop/PRACTICAS/data_bladder")

# ---------------------------
# 2. Load final Seurat object and external metadata
# ---------------------------
seurat <- readRDS(paste0(sample_id, "_final.rds"))
metadata_full <- read.delim("metadata_full.txt", stringsAsFactors = FALSE)

# ---------------------------
# 3. Filter metadata for cancer spots with clone info
# ---------------------------
metadata_filtered <- metadata_full %>%
  filter(sample == sample_id,
         new_compartment == "Cancer",
         !is.na(clone_grouped))

# ---------------------------
# 4. Sanity check: barcode match
# ---------------------------
barcodes <- metadata_filtered$spotid
barcodes_seurat <- intersect(barcodes, colnames(seurat))
if (length(barcodes_seurat) == 0) {
  stop("No matching barcodes found in Seurat object!")
} else {
  cat("Found", length(barcodes_seurat), "matching barcodes.\n")
}

# ---------------------------
# 5. Subset Seurat object to those barcodes
# ---------------------------
seurat_filtered <- subset(seurat, cells = barcodes_seurat)

# ---------------------------
# 6. Set RNA as default assay
# ---------------------------
DefaultAssay(seurat_filtered) <- "RNA"
if (DefaultAssay(seurat_filtered) != "RNA") {
  stop("Default assay is not RNA!")
} else {
  cat("Default assay is RNA.\n")
}

# ---------------------------
# 7. Define Y chromosome genes
# ---------------------------
y_genes <- c("AMELY", "BPY2", "BPY2B", "BPY2C", "CDY1A", "CDY1B", "CDY2A", "CDY2B",
             "DAZ1", "DAZ2", "DAZ3", "DAZ4", "DDX3Y", "EIF1AY", "GOLGA2LY", "HSFY1",
             "HSFY2", "KDM5D", "PCDH11Y", "PRY", "RBMY1A1", "RBMY1B", "RBMY1C",
             "RBMY1D", "RBMY1E", "RBMY1F", "RBMY1J", "RPS4Y1", "RPS4Y2", "SRY",
             "TBL1Y", "TGIF2LY", "TMSB4Y", "TXLNGY", "TTY1", "TTY2", "TTY3", "TTY4",
             "USP9Y", "UTY", "VCY", "XKRY", "ZFY", "BCORP1")

# ---------------------------
# 8. Get detected Y genes and calculate raw sums
# ---------------------------
y_genes_detected <- intersect(y_genes, rownames(seurat_filtered))
cat("Detected Y genes:", paste(y_genes_detected, collapse = ", "), "\n")

# Compute Y_sum (raw counts sum per spot)
Y_sum <- colSums(as.matrix(seurat_filtered@assays$RNA@counts[y_genes_detected, , drop = FALSE]), na.rm = TRUE)

# ---------------------------
# 9. Join Y_sum with metadata
# ---------------------------
metadata_plot <- metadata_filtered %>%
  filter(spotid %in% names(Y_sum)) %>%
  mutate(Y_sum = Y_sum[spotid])


# ---------------------------
# 10. Boxplot of ChrY raw counts per clone + threshold line
# ---------------------------
ggplot(metadata_plot, aes(x = clone_grouped, y = Y_sum, fill = clone_grouped)) +
  geom_boxplot(alpha = 0.7, outlier.color = "black") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red", size = 1) +  # ← Threshold here
  labs(title = paste("Raw ChrY Expression per Clone -", sample_id),
       subtitle = "Red dashed line = spot-level ChrY loss threshold (Y_sum < 5)",
       x = "Clone Group", y = "Y_sum (ChrY Total Counts)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---------------------------
# 12. Add binary loss call per spot
# ---------------------------
metadata_plot$ChrY_loss <- metadata_plot$Y_sum < 5

# ---------------------------
# 13. Summarize ChrY loss % per clone
# ---------------------------
clone_summary <- metadata_plot %>%
  group_by(clone_grouped) %>%
  summarise(
    total_spots = n(),
    chrY_loss_spots = sum(ChrY_loss),
    chrY_loss_percent = round((chrY_loss_spots / total_spots) * 100, 1)
  ) %>%
  mutate(ChrY_loss_clone_level = ifelse(chrY_loss_percent >= 30, "ChrY_LOSS", "ChrY_PRESENT"))

# Optional: print or exportation
print(clone_summary)
write.csv(clone_summary, paste0(sample_id, "_ChrY_loss_by_clone.csv"), row.names = FALSE)


# ---------------------------
# 14. Export result table
# ---------------------------
write.csv(metadata_plot, file = paste0(sample_id, "_Ysum_clones.csv"), row.names = FALSE)
