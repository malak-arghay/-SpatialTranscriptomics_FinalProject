# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# ---------------------------
# Set working directory and sample
# ---------------------------
setwd("~/Desktop/PRACTICAS/data_bladder")
sample_id <- "DU34"

# ---------------------------
# Load Seurat object and metadata
# ---------------------------
data_single <- readRDS(paste0(sample_id, ".rds"))
data_single@images[[1]]@scale.factors$spot <- 80
metadata_full <- read.delim("metadata_full.txt", stringsAsFactors = FALSE)

# ---------------------------
# Subset to TME spots only
# ---------------------------
barcodes_tme <- metadata_full$spotid[metadata_full$sample == sample_id &
                                       metadata_full$new_compartment == "TME"]
data_tme <- subset(data_single, cells = barcodes_tme)
cat("TME spots in sample", sample_id, ":", ncol(data_tme), "\n")

# ---------------------------
# Define Y chromosome genes
# ---------------------------
y_genes <- c("AMELY", "BPY2", "BPY2B", "BPY2C", "CDY1A", "CDY1B", "CDY2A", "CDY2B",
             "DAZ1", "DAZ2", "DAZ3", "DAZ4", "DDX3Y", "EIF1AY", "GOLGA2LY", "HSFY1",
             "HSFY2", "KDM5D", "PCDH11Y", "PRY", "RBMY1A1", "RBMY1B", "RBMY1C",
             "RBMY1D", "RBMY1E", "RBMY1F", "RBMY1J", "RPS4Y1", "RPS4Y2", "SRY",
             "TBL1Y", "TGIF2LY", "TMSB4Y", "TXLNGY", "TTY1", "TTY2", "TTY3", "TTY4",
             "USP9Y", "UTY", "VCY", "XKRY", "ZFY", "BCORP1")

# ---------------------------
# Calculate raw sum of ChrY expression
# ---------------------------
DefaultAssay(data_tme) <- "RNA"
y_genes_detected <- intersect(y_genes, rownames(data_tme))
cat(" genes detected in TME:", paste(y_genes_detected, collapse = ", "), "\n")

if (length(y_genes_detected) > 0) {
  counts_matrix <- data_tme@assays$RNA@counts
  y_sum <- colSums(as.matrix(counts_matrix[y_genes_detected, , drop = FALSE]), na.rm = TRUE)
  data_tme$Y_sum <- y_sum
  data_tme$ChrY_loss <- y_sum < 5 
  
  # ---------------------------
  # Spatial plot: raw ChrY sum in TME
  # ---------------------------
  SpatialFeaturePlot(data_tme, features = "Y_sum") +
    ggtitle(paste("Raw ChrY Sum in TME -", sample_id)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # ---------------------------
  # Spatial plot: ChrY loss (binary)
  # ---------------------------
  data_tme$ChrY_loss_label <- ifelse(data_tme$ChrY_loss, "Loss", "Expression")
  Idents(data_tme) <- "ChrY_loss_label"
  
  SpatialDimPlot(data_tme, label = TRUE, repel = TRUE,
                 cols = c("Loss" = "black", "Expression" = "skyblue")) +
    ggtitle(paste("ChrY Loss Map in TME -", sample_id)) +
    theme(plot.title = element_text(hjust = 0.5))
  
} else {
  cat("No Y genes detected in this TME subset.\n")
}

total_spots <- ncol(data_tme)
loss_spots <- sum(data_tme$ChrY_loss)
percent_loss <- round((loss_spots / total_spots) * 100, 1)
cat("ChrY loss in TME (", sample_id, "):", percent_loss, "% of spots\n")

