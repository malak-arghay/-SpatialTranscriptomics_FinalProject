## 1. Are clones with ChrY_LOSS more likely to be non-responders (NR)? 

## merging reponse into merged_all
any(duplicated(clone_classes[, c("sample", "clone")]))

merged_all <- merged_all %>%
  left_join(
    clone_classes %>%
      select(sample, clone_grouped = clone, response_combined),
    by = c("sample", "clone_grouped")
  )

table(is.na(merged_all$response_combined))

# STEP 1: Build a contingency table

# Remove NAs just to be clean
table_chrY_response <- merged_all %>%
  filter(!is.na(response_combined), !is.na(ChrY_loss_call)) %>%
  count(ChrY_loss_call, response_combined) %>%
  tidyr::pivot_wider(names_from = response_combined, values_from = n, values_fill = 0)

print(table_chrY_response)

#STEP 2: Statistical test

# Re-create contingency table for test
contingency_test <- table(merged_all$ChrY_loss_call, merged_all$response_combined)

# Run the test
chisq.test(contingency_test)

## STEP 3: Optional — Bar plot (clean visual)
library(dplyr)
library(ggplot2)

# Prepare data
merged_all %>%
  filter(!is.na(response_combined)) %>%
  group_by(ChrY_loss_call, response_combined) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ChrY_loss_call) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(aes(x = ChrY_loss_call, y = freq, fill = response_combined)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("CR/PR" = "#1b9e77", "NR" = "#d95f02", "Missing" = "grey")) +
  labs(
    title = "Treatment Response Distribution by ChrY Loss Status",
    x = "ChrY Loss Status", y = "Proportion of Clones",
    fill = "Response"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11))



##2.  Loop for all subtype-specific
results <- list()

# Get unique subtypes (removing NA just in case)
subtypes <- unique(na.omit(merged_all$consensusClass))

for (subtype in subtypes) {
  # Filter to current subtype
  subset_data <- merged_all %>%
    filter(consensusClass == subtype)
  
  # Create contingency table: ChrY loss vs. response_combined
  tab <- table(subset_data$ChrY_loss_call, subset_data$response_combined)
  
  # Only run test if table has >= 2 rows and columns
  if (nrow(tab) >= 2 && ncol(tab) >= 2) {
    test_result <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
    
    # Store results
    results[[subtype]] <- list(
      table = tab,
      p_value = test_result$p.value
    )
  } else {
    # Store reason for skipping
    results[[subtype]] <- list(
      table = tab,
      p_value = NA,
      note = "Skipped: Not enough categories to run test"
    )
  }
}

# Print nicely
for (subtype in names(results)) {
  cat("\nSubtype:", subtype, "\n")
  print(results[[subtype]]$table)
  cat("P-value:", results[[subtype]]$p_value, "\n")
  if (!is.null(results[[subtype]]$note)) {
    cat("Note:", results[[subtype]]$note, "\n")
  }
}


# ───────────────────────────────────────────────────────────────
# FINAL VISUALIZATION: Faceted bar plot of response by ChrY loss per subtype
# ───────────────────────────────────────────────────────────────

# 1. Prepare data: remove NA or 'Missing' responses and keep only desired subtypes
plot_data <- merged_all %>%
  filter(
    !is.na(consensusClass),
    !is.na(ChrY_loss_call),
    !is.na(response_combined),
    response_combined != "Missing",
    consensusClass %in% c("LumP", "LumNS", "LumU")
  )

# 2. Create faceted bar plot showing response distribution per ChrY status within subtypes
ggplot(plot_data, aes(x = ChrY_loss_call, fill = response_combined)) +
  geom_bar(position = "fill") +
  facet_wrap(~ consensusClass) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = c("CR/PR" = "#1b9e77", "NR" = "#d95f02")
  ) +
  labs(
    title = "Treatment Response by ChrY Loss Status within Molecular Subtypes",
    subtitle = "Proportion of clones with CR/PR vs NR by ChrY status in LumP, LumNS, and LumU",
    x = "ChrY Loss Status",
    y = "Proportion of Clones",
    fill = "Response"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




            