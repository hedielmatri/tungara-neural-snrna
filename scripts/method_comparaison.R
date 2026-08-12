library(qs)
library(pheatmap)
library(dplyr)
library(mclust)
setwd("/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/")

frog_jace = qread("Jace_R/combined_int_test.qs")
frog_hedi = readRDS("Hedi_R/tunga_frog_blast_with_gaba_subtypes.rds")

meta_jace = data.frame(
  cell_id = colnames(frog_jace), # rownames(frog_jace@meta.data)
  label_jace = frog_jace$sctype_label_45PC_0.6
)

meta_hedi = data.frame(
  cell_id = colnames(frog_hedi), # rownames(frog_hedi@meta.data)
  label_hedi = frog_hedi$refined_annotation
)

merged_data = merge(meta_jace, meta_hedi, by = "cell_id")

comparison_table = table(merged_data$label_jace, merged_data$label_hedi)
#print(comparison_table)


prob_table = prop.table(comparison_table, margin = 1) 
pheatmap(prob_table,  
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = TRUE, 
         main = "Jace against Hedi")


head(colnames(frog_jace))
head(colnames(frog_hedi))

shared_cells = intersect(colnames(frog_jace), colnames(frog_hedi))
length(shared_cells)


rowSums(comparison_table)


pheatmap(comparison_table, display_numbers = TRUE, main = "Raw Counts")



best_matches = apply(comparison_table, 1, max)
agreement = sum(best_matches) / sum(comparison_table)
print(paste0("Agreement:", round(agreement * 100, 2), "%"))


meta_jace = data.frame(
  cell_id = colnames(frog_jace), # rownames(frog_jace@meta.data)
  label_jace = frog_jace$sctype_label_45PC_0.6
)
meta_hedi = data.frame(
  cell_id = colnames(frog_hedi), # rownames(frog_hedi@meta.data)
  label_hedi = frog_hedi$refined_annotation
)
merged_data = merge(meta_jace, meta_hedi, by = "cell_id")
comparison_table = table(merged_data$label_jace, merged_data$label_hedi)

best_matches = apply(comparison_table, 1, max)
agreement = sum(best_matches) / sum(comparison_table)
print(paste0("Agreement:", round(agreement * 100, 2), "%"))



ari_score = adjustedRandIndex(merged_data$label_jace, merged_data$label_hedi)
print(paste0("Adjusted rand index: ", round(ari_score, 3)))


# install.packages("ggalluvial")
library(ggalluvial)
library(ggplot2)

plot_data = as.data.frame(table(merged_data$label_jace, merged_data$label_hedi))
colnames(plot_data) = c("Jace", "Hedi", "Freq")

# Plot the Sankey / Alluvial diagram
ggplot(plot_data, aes(axis1 = Jace, axis2 = Hedi, y = Freq)) +
  geom_alluvium(aes(fill = Jace), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Jace's Labels", "Hedi's Labels")) +
  theme_minimal() +
  ggtitle("Cell Label Mapping: Jace vs. Hedi") +
  theme(legend.position = "none")





merged_data_broad = merged_data %>%
  mutate(
    jace_broad = case_when(
      grepl("Excitatory neuron", label_jace, ignore.case = TRUE) ~ "Excitatory Neuron",
      grepl("GABAergic neuron", label_jace, ignore.case = TRUE) ~ "GABAergic Neuron",
      grepl("astrocyte", label_jace, ignore.case = TRUE) ~ "Astrocyte",
      grepl("oligodendro", label_jace, ignore.case = TRUE) ~ "Oligodendrocyte",
      grepl("endothelial", label_jace, ignore.case = TRUE) ~ "Endothelial Cell",
      grepl("antigen-presenting|microglia|immune", label_jace, ignore.case = TRUE) ~ "Immune Cell",
      TRUE ~ label_jace
    ),
    hedi_broad = case_when(
      grepl("Excitatory neuron", label_hedi, ignore.case = TRUE) ~ "Excitatory Neuron",
      grepl("GABAergic neuron", label_hedi, ignore.case = TRUE) ~ "GABAergic Neuron",
      grepl("astrocyte", label_hedi, ignore.case = TRUE) ~ "Astrocyte",
      grepl("oligodendro", label_hedi, ignore.case = TRUE) ~ "Oligodendrocyte",
      grepl("endothelial", label_hedi, ignore.case = TRUE) ~ "Endothelial Cell",
      grepl("antigen-presenting|microglia|immune", label_hedi, ignore.case = TRUE) ~ "Immune Cell",
      TRUE ~ label_hedi
    )
  )

# Filter for ONLY DISAGREEments
disagreements_only <- merged_data_broad %>%
  filter(jace_broad != hedi_broad)

plot_data = as.data.frame(table(disagreements_only$jace_broad, disagreements_only$hedi_broad))
colnames(plot_data) = c("Jace", "Hedi", "Freq")

# Remove 0 frequency rows to keep the plot clean
plot_data = plot_data %>% filter(Freq > 0)

ggplot(plot_data, aes(axis1 = Jace, axis2 = Hedi, y = Freq)) +
  geom_alluvium(aes(fill = Jace), width = 1/12, alpha = 0.7) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  scale_x_discrete(limits = c("Jace's labels", "Hedi's labels")) +
  theme_minimal() +
  labs(
    title = "Where the pipelines diverged",
    subtitle = paste("Showing the", nrow(disagreements_only), "cells with conflicting broad annotations"),
    y = "Number of cells"
  ) +
  theme(legend.position = "none")


comparison_table = table(merged_data_broad$label_jace, merged_data_broad$label_hedi)
best_matches = apply(comparison_table, 1, max)
agreement = sum(best_matches) / sum(comparison_table)
print(paste0("Agreement:", round(agreement * 100, 2), "%"))



ari_score = adjustedRandIndex(merged_data_broad$label_jace, merged_data_broad$label_hedi)
print(paste0("Adjusted rand index: ", round(ari_score, 3)))

 