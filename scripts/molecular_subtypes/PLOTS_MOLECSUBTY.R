# 1st plot
# Average raw ChrY expression per clone (avg_Ysum), grouped by molecular subtype
library(ggplot2)

ggplot(merged_all %>% filter(!is.na(consensusClass)), aes(x = consensusClass, y = avg_Ysum, fill = consensusClass)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "ChrY Expression (avg_Ysum) per Molecular Subtype",
       x = "Molecular Subtype", y = "Average ChrY Expression per Clone") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##2nd plot , spot level
ggplot(merged_all %>% filter(!is.na(consensusClass)), aes(x = consensusClass, y = percent_loss, fill = consensusClass)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "ChrY Loss % per Molecular Subtype",
       x = "Molecular Subtype", y = "ChrY Loss (% of spots)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## clone level

library(ggplot2)
library(dplyr)

# Filter NAs just in case
merged_all_clean <- merged_all %>% filter(!is.na(consensusClass), !is.na(ChrY_loss_call))

# Count clones per subtype + loss call
clone_counts <- merged_all_clean %>%
  group_by(consensusClass, ChrY_loss_call) %>%
  summarise(n_clones = n(), .groups = "drop")

# Bar plot of number of clones per subtype, by ChrY loss status
ggplot(clone_counts, aes(x = consensusClass, y = n_clones, fill = ChrY_loss_call)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("ChrY_LOSS" = "red", "ChrY_PRESENT" = "grey")) +
  labs(
    title = "Number of Clones per Molecular Subtype by ChrY Loss Status",
    x = "Molecular Subtype", y = "Number of Clones", fill = "ChrY Loss Status"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##3rd plot
library(dplyr)

merged_all %>%
  group_by(consensusClass, ChrY_loss_call) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(aes(x = consensusClass, y = freq, fill = ChrY_loss_call)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "ChrY Loss Call Distribution per Subtype",
       x = "Molecular Subtype", y = "Proportion",
       fill = "ChrY Loss Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



##Staistics for molecular subtypes

##1st statistical test
#Kruskal–Wallis Test del reves 

library(dplyr)

# Prepare data (ensure no NAs)
kw_data <- merged_all %>%
  filter(!is.na(consensusClass)) %>%
  select(consensusClass, percent_loss)

# Kruskal–Wallis
kw_result <- kruskal.test(percent_loss ~ consensusClass, data = kw_data)
print(kw_result)
##highly significant


pairwise_results <- pairwise.wilcox.test(
  x = kw_data$percent_loss,
  g = kw_data$consensusClass,
  p.adjust.method = "BH"
)
print(pairwise_results)


##2 statisticsal test 
#Chi-square test of independence 

##1 . Build the contingency table
chrY_contingency <- table(merged_all$ChrY_loss_call, merged_all$consensusClass)
# View the table
chrY_contingency

##2.Run the Chi-square test
chisq_test <- chisq.test(chrY_contingency)

# Print results
print(chisq_test)



























