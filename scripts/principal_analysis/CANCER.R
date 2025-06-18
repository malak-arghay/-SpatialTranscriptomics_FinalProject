# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# ---------------------------
# Setting working directory and sample
# ---------------------------
setwd("~/Desktop/PRACTICAS/data_bladder")
sample_id <- "DU2"  # ← change for each sample

# ---------------------------
# Load Seurat object and metadata
# ---------------------------
data_single <- readRDS(paste0(sample_id, ".rds"))
data_single@images[[1]]@scale.factors$spot <- 80
metadata_full <- read.delim("metadata_full.txt", stringsAsFactors = FALSE)

# ---------------------------
# Sanity check barcodes
# ---------------------------
meta_sub <- metadata_full[metadata_full$sample == sample_id, c("spotid", "new_compartment")]
rownames(meta_sub) <- meta_sub$spotid
common_barcodes <- intersect(colnames(data_single), meta_sub$spotid)
if (length(common_barcodes) == 0) stop("No overlapping barcodes found.") else
  cat(length(common_barcodes), "matching barcodes found.\n")

# ---------------------------
# Subset Seurat object to all barcodes
# ---------------------------
barcodes_all <- meta_sub$spotid
data_all_spots <- subset(data_single, cells = barcodes_all)
meta_sub <- meta_sub[colnames(data_all_spots), ]
data_all_spots$compartment_definitive <- meta_sub$new_compartment

# ---------------------------
# Plot compartments
# ---------------------------
SpatialDimPlot(data_all_spots, group.by = "compartment_definitive",
               cols = c("Cancer" = "red", "TME" = "green", "Buffer" = "lightblue")) +
  ggtitle(paste("Compartments Overview -", sample_id)) +
  theme(plot.title = element_text(hjust = 0.5))

# ---------------------------
# Subset to cancer spots only
# ---------------------------
barcodes_cancer <- metadata_full$spotid[metadata_full$sample == sample_id &
                                          metadata_full$new_compartment == "Cancer"]
data_cancer <- subset(data_single, cells = barcodes_cancer)
cat("Cancer spots in sample", sample_id, ":", ncol(data_cancer), "\n")

# View cancer layout
Idents(data_cancer) <- "Cancer"
SpatialFeaturePlot(data_cancer, features = NULL) +
  ggtitle(paste("Tissue Plot -", sample_id))


# ---------------------------
# Define Y genes (HGNC-approved list)
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
DefaultAssay(data_cancer) <- "RNA"
y_genes_detected <- intersect(y_genes, rownames(data_cancer))
cat("Y genes detected:", paste(y_genes_detected, collapse = ", "), "\n")

if (length(y_genes_detected) > 0) {
  counts_matrix <- data_cancer@assays$RNA@counts
  y_sum <- colSums(as.matrix(counts_matrix[y_genes_detected, , drop = FALSE]), na.rm = TRUE)
  data_cancer$Y_sum <- y_sum
  data_cancer$ChrY_loss <- y_sum < 5
  # ---------------------------
  # Summary: ChrY loss % in cancer spots
  # ---------------------------
  total_cancer_spots <- ncol(data_cancer)
  loss_cancer_spots <- sum(data_cancer$ChrY_loss)
  percent_cancer_loss <- round((loss_cancer_spots / total_cancer_spots) * 100, 1)
  
  cat("ChrY loss % in cancer spots:", percent_cancer_loss, "%\n")
  
  # ---------------------------
  #  Spatial plot of raw ChrY sum
  # ---------------------------
  SpatialFeaturePlot(data_cancer, features = "Y_sum") +
    ggtitle(paste("Raw Sum of ChrY Genes -", sample_id))
  
  # ---------------------------
  #  Spatial plot of ChrY loss (binary)
  # ---------------------------
  data_cancer$ChrY_loss_label <- ifelse(data_cancer$ChrY_loss, "Loss", "Expression")
  Idents(data_cancer) <- "ChrY_loss_label"
  
  SpatialDimPlot(data_cancer, label = TRUE, repel = TRUE,
                 cols = c("Loss" = "black", "Expression" = "skyblue")) +
    ggtitle(paste("ChrY Loss Map -", sample_id)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # ---------------------------
  # Filter metadata for clone analysis
  # ---------------------------
  metadata_filtered <- metadata_full %>%
    filter(sample == sample_id,
           new_compartment == "Cancer",
           !is.na(clone_grouped))
  
  # Attach Y_sum and ChrY_loss to metadata
  metadata_filtered$Y_sum <- y_sum[metadata_filtered$spotid]
  metadata_filtered$ChrY_loss <- metadata_filtered$Y_sum < 5
  
  # ---------------------------
  # Count loss per clone
  # ---------------------------
  clone_summary <- metadata_filtered %>%
    group_by(clone_grouped) %>%
    summarise(
      total_spots = n(),
      chrY_loss_spots = sum(ChrY_loss),
      chrY_loss_percent = round((chrY_loss_spots / total_spots) * 100, 1)
    )
  
  # ---------------------------
  # Add binary call per clone based on threshold
  # ---------------------------
  threshold <- 30  # Threshold in percentage (you can change it later if needed)
  clone_summary$ChrY_loss_clone_level <- ifelse(clone_summary$chrY_loss_percent >= threshold,
                                                "ChrY_LOSS", "ChrY_PRESENT")
  
  
  # Plot fraction of ChrY loss per clone
  ggplot(clone_summary, aes(x = clone_grouped, y = chrY_loss_percent, fill = clone_grouped)) +
    geom_col() +
    labs(title = paste("ChrY Loss % per Clone -", sample_id),
         x = "Clone", y = "ChrY Loss (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #Plot threashold
  ggplot(clone_summary, aes(x = clone_grouped, y = chrY_loss_percent,
                            fill = ChrY_loss_clone_level)) +
    geom_col() +
    scale_fill_manual(values = c("ChrY_LOSS" = "red", "ChrY_PRESENT" = "grey")) +
    labs(title = paste("ChrY Loss % per Clone -", sample_id),
         x = "Clone", y = "ChrY Loss (%)", fill = "ChrY Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  # Save output
  write.csv(clone_summary, paste0(sample_id, "_ChrY_loss_by_clone.csv"), row.names = FALSE)
  
} else {
  cat("❌ No ChrY genes detected in", sample_id, "\n")
}

