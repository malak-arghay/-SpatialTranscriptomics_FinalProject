## CODE FOR INDIVIDUAL ANALSIYS OF EACH MOLEC SUBTYPE + Y SCORE BINNARY

# ---------------------------
# Set sample and working directory
# ---------------------------
sample_id <- "DU**"
setwd("~/Desktop/PRACTICAS/data_bladder")

# ---------------------------
# Load metadata with Y_sum and ChrY_loss per spot
# ---------------------------
metadata_plot <- read.csv(paste0(sample_id, "_Ysum_clones.csv"))

# ---------------------------
# Load molecular subtype file
# ---------------------------
clone_classes <- read.delim("clone_classes.txt", stringsAsFactors = FALSE)

# ---------------------------
# Filter molecular subtype info for the current sample
# ---------------------------
subtype_sample <- clone_classes %>%
  filter(sample == sample_id)

# ---------------------------
# Compute per-clone summary from the metadata
# ---------------------------
clone_summary <- metadata_plot %>%
  group_by(clone_grouped) %>%
  summarise(
    total_spots = n(),
    avg_Ysum = round(mean(Y_sum, na.rm = TRUE), 2),
    percent_loss = round((sum(Y_sum < 5) / n()) * 100, 1)
  ) %>%
  mutate(
    ChrY_loss_call = ifelse(percent_loss >= 30, "ChrY_LOSS", "ChrY_PRESENT")
  )

# 2. Load and filter molecular subtypes
subtype_sample <- clone_classes %>% filter(sample == sample_id)

# Sanity check before merging
cat("\n Clones in Y-summary file (Seurat):\n", 
    paste(sort(unique(clone_summary$clone_grouped)), collapse = ", "), "\n\n")

cat(" Clones in molecular subtype file (", sample_id, "):\n", 
    paste(sort(unique(subtype_sample$clone)), collapse = ", "), "\n\n")

# ---------------------------
# Join with molecular subtype annotation
# ---------------------------
clone_summary <- clone_summary %>%
  left_join(subtype_sample, by = c("clone_grouped" = "clone"))

# ---------------------------
# Save output
# ---------------------------
write.csv(clone_summary, paste0(sample_id, "_ChrYloss_vs_Subtype.csv"), row.names = FALSE)

# Done
cat("Summary table created for", sample_id, "\n")

