# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# ---------------------------
# 1. Sample ID and working directory
# ---------------------------
sample_id <- "DU34"
setwd("~/Desktop/PRACTICAS/data_bladder")

# ---------------------------
# 2. Load Seurat object and metadata
# ---------------------------
seurat <- readRDS(paste0(sample_id, "_final.rds"))
seurat@images[[1]]@scale.factors$spot <- 80  # adjust as needed for display

metadata_full <- read.delim("metadata_full.txt", stringsAsFactors = FALSE)

# ---------------------------
# 3. Filter external metadata for this sample and cancer spots with known clone
# ---------------------------
metadata_filtered <- metadata_full %>%
  filter(sample == sample_id,
         new_compartment == "Cancer",
         !is.na(clone_grouped))

# ---------------------------
# 4. Get valid barcodes that exist in both metadata and Seurat object
# ---------------------------
barcodes <- intersect(metadata_filtered$spotid, colnames(seurat))
if (length(barcodes) == 0) {
  stop("No matching barcodes found between Seurat object and metadata.")
} else {
  cat("Sanity check passed:", length(barcodes), "matching barcodes found.\n")
}

# ---------------------------
# 5. Subset Seurat object to valid cancer+clone barcodes
# ---------------------------
seurat_filtered <- subset(seurat, cells = barcodes)

# ---------------------------
# 6. Ensure barcode order matches between Seurat and metadata
# ---------------------------
metadata_filtered <- metadata_filtered[match(colnames(seurat_filtered), metadata_filtered$spotid), ]

if (!all(metadata_filtered$spotid == colnames(seurat_filtered))) {
  stop("Barcode order mismatch between metadata and Seurat object!")
} else {
  cat("Barcode order matches.\n")
}

# ---------------------------
# 7. Add clone_grouped info to Seurat metadata (from trusted external file)
# ---------------------------
seurat_filtered$clone_grouped <- metadata_filtered$clone_grouped
Idents(seurat_filtered) <- seurat_filtered$clone_grouped

# ---------------------------
# 8. Plot spatial distribution of clones in the tissue
# ---------------------------
SpatialDimPlot(seurat_filtered, group.by = "clone_grouped", label = TRUE, repel = TRUE, label.size = 5) +
  ggtitle(paste("Spatial Clone Map -", sample_id)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(sample_id, "_clone_spatial_map.png"), width = 6, height = 5)




