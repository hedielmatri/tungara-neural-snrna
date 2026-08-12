setwd("/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/Hedi_R")

library(Seurat)
library(dplyr)
library(readxl)
library(tidyr)
library(HGNChelper)
library(ggplot2)
library(future)
library(devtools)
library(hdWGCNA)
library(WGCNA)
library(writexl)
library(pheatmap)

options(future.globals.maxSize = 64000 * 1024^2)
Sys.setenv(R_MAX_VSIZE = "64Gb")





tunga_frog_merged_blast = readRDS("tunga_frog_annotated_blast_wgcna.rds")
neural_obj = readRDS("neural_wgcna_obj.rds")
translation_table = read_excel("translation_table_xenopus_to_tungara.xlsx")


all_response_genes = c("NRTN", "CAPRIN2", "NUMB", "RPGRIP1L", "SLC22A2", "ATP7B", "TRAF2", "SP4", "SMAD2", "E2F3", "MAPK8", "AVP")




# FIND WHAT CLUSTERS THOSE GENES ARE MOST EXPRESSED IN


# Adding module score to every cell based on 12 response genes
neural_obj = AddModuleScore(
  object = neural_obj,
  features = list(all_response_genes),
  name = "WhineChuck_response_score"
)

FeaturePlot(
  neural_obj, features = "WhineChuck_response_score1", split.by = "orig.ident", reduction = "umap_by_pca",
            order = TRUE,
            pt.size = .5) + 
  labs(title = "Localization of all response genes")


score_data = FetchData(
  object = neural_obj, 
  vars = c("WhineChuck_response_score1", "refined_annotation", "orig.ident")
)

aggregate_table = score_data %>%
  group_by(refined_annotation, orig.ident) %>%
  summarize(Mean_response_score = mean(WhineChuck_response_score1, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Mean_response_score))

print(aggregate_table)







# DEG ANALYSIS ON ONLY THOSE CELLS


substantive_cells = subset(neural_obj, subset = refined_annotation %in% c("GABAergic neuron_tac1 high", 
                                                                         "Excitatory neuron_tac1 high", 
                                                                         "GABAergic neuron_sst high"))
Idents(substantive_cells) = "orig.ident"

# DEG Whine-chuck against Whine
diff_expr_genes = FindMarkers(substantive_cells, 
                                     ident.1 = "Whine-chuck", 
                                     ident.2 = "Whine", 
                                     logfc.threshold = 0.25, 
                                     min.pct = 0.1)

significant_upregulated = diff_expr_genes %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

print(head(significant_upregulated, 100))
print(head(translation_table))


# TRANSLATE SOME LOC GENES

significant_upregulated$tungara_symbol = rownames(significant_upregulated)

translated_results = left_join(significant_upregulated, translation_table, by = "tungara_symbol")

final_table = translated_results %>%
  select(tungara_symbol, xenopus_symbol, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  arrange(desc(avg_log2FC))

print(head(final_table, 30))





# CHECK Feature DIVIDe BETWEEN CApRIN2 AND NRTN

# looking at Whine-chuck substantive cells
wc_substantive = subset(substantive_cells, subset = orig.ident == "Whine-chuck")

FeatureScatter(wc_substantive, 
               feature1 = "CAPRIN2", 
               feature2 = "NRTN", 
               group.by = "refined_annotation", 
               pt.size = 1.5) +
  labs(title = "Endocrine against structural ",
       x = "CAPRIN2 Expression: Hormone",
       y = "NRTN Expression: Rewiring")



expression_data = FetchData(wc_substantive, vars = c("CAPRIN2", "NRTN"))

expression_data = expression_data %>%
  mutate(
    strategy = case_when(
      CAPRIN2 > 0 & NRTN > 0 ~ "General_both",
      CAPRIN2 > 0 & NRTN == 0 ~ "endocrine_only",
      CAPRIN2 == 0 & NRTN > 0 ~ "structural_only",
      TRUE ~ "Inactive"
    )
  )

wc_substantive = AddMetaData(wc_substantive, metadata = expression_data$strategy, col.name = "strategy")

# split sizes
table(wc_substantive$strategy)



Idents(wc_substantive) = "strategy"

# Finding drivers of endocrine and structural
strategy_markers = FindMarkers(wc_substantive, 
                                ident.1 = "endocrine_only", 
                                ident.2 = "structural_only", 
                                logfc.threshold = 0.25)

strategy_markers$tungara_symbol = rownames(strategy_markers)
translated_results2 = left_join(strategy_markers, translation_table, by = "tungara_symbol")

strategy_table = translated_results2 %>%
  select(tungara_symbol, xenopus_symbol, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

# Positive log2fc = endocrine
# Negative log2fc = structural
print(head(strategy_table, 20))
print(tail(strategy_table, 20))



##### I ADDED INDEPENDENTLY GENOTYPE INFORMATION!
# !!!!!

genotype_strategy = FetchData(wc_substantive, vars = c("strategy", "souporcell_assignment"))

# proportion of each strategy used by each frog
genotype_strat = genotype_strategy %>%
  filter(strategy != "Inactive") %>% 
  group_by(souporcell_assignment, strategy) %>%
  summarize(Cell_count = n(), .groups = "drop") %>%
  group_by(souporcell_assignment) %>%
  mutate(Strategy_percent = round((Cell_count / sum(Cell_count)) * 100, 2)) %>%
  arrange(souporcell_assignment, desc(Strategy_percent))

# Does genetics control the brain's strategy?
print(as.data.frame(genotype_strat))


###########

# Finding correlated genes

expression_matrix = as.matrix(GetAssayData(wc_substantive, assay = "RNA", layer = "data"))

caprin2_cor = cor(expression_matrix["CAPRIN2", ], t(expression_matrix))
caprin2_top = as.data.frame(sort(caprin2_cor[1, ], decreasing = TRUE)[2:15]) # Top 14 correlated
colnames(caprin2_top) = "corr" 

nrtn_cor = cor(expression_matrix["NRTN", ], t(expression_matrix))
nrtn_top = as.data.frame(sort(nrtn_cor[1, ], decreasing = TRUE)[2:15])
colnames(nrtn_top) = "corr" 


caprin2_top$tungara_symbol = rownames(caprin2_top)
translated_results3 = left_join(caprin2_top, translation_table, by = "tungara_symbol")

caprin2_top = translated_results3 %>%
  select(tungara_symbol, xenopus_symbol, corr)
print(caprin2_top)


nrtn_top$tungara_symbol = rownames(nrtn_top)
translated_results4 = left_join(nrtn_top, translation_table, by = "tungara_symbol")

nrtn_top = translated_results4 %>%
  select(tungara_symbol, xenopus_symbol, corr)
print(nrtn_top)


caprin2_besthit = c("CAPRIN2", "DNASE1", "JARID2", "DMPK", "SLC4A7", "DEDD2", "PNISR", "TRAF2", "GRIN2D", "NVL")

VlnPlot(wc_substantive, features = caprin2_besthit, group.by = "strategy", pt.size = 0.1, ncol = 5)


nrtn_besthit = c("NRTN", "SHISA8", "DDAH1", "WIF1", "ATP7B", "HIP1R", "PRKDC", "NUMB", "LOC140071013", "HMGCLL1")

VlnPlot(wc_substantive, features = nrtn_besthit, group.by = "strategy", pt.size = 0.1, ncol = 5) 
# "LOC140071013 = LMNTD2"






# FIND RECEPTORS AND MAKERS

wc_substantive$Response_state = ifelse(
  wc_substantive$strategy == "Inactive", 
  "Inactive", 
  "Active"
)

Idents(wc_substantive) = "Response_state"

# Finding genes driving active state
active_genes = FindMarkers(wc_substantive, 
                             ident.1 = "Active", 
                             ident.2 = "Inactive", 
                             logfc.threshold = 0.25)


sgnificant_active_genes = active_genes %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

sgnificant_active_genes$tungara_symbol = rownames(sgnificant_active_genes)
translated_results5 = left_join(sgnificant_active_genes, translation_table, by = "tungara_symbol")

sgnificant_active_genes = translated_results5 %>%
  select(tungara_symbol, xenopus_symbol, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  head(30) %>%
  print()

# top 30 differences between a cell that fires and a cell that is inactive



# FIND RECEPTORS AND MAKERS SECOND TRIAL

baseline_substantive_cells = subset(substantive_cells, subset = orig.ident == "Whine")

baseline_substantive_expression = GetAssayData(baseline_substantive_cells, assay = "RNA", layer = "data")

# Calculating what percentage of cells express each gene at baseline
baseline_substantive_pct = rowSums(baseline_substantive_expression > 0) / ncol(baseline_substantive_expression)

# Filtering high baseline expression expressed in >50% of the cells
primed_substantive_genes = sort(baseline_substantive_pct[baseline_substantive_pct > 0.5], decreasing = TRUE)

primed_substantive_df = data.frame(gene = names(primed_substantive_genes), baseline_substantive_pct = primed_substantive_genes)
print(head(primed_substantive_df, 20))
print(sum(primed_substantive_df$baseline_substantive_pct > 0.5)/nrow(primed_substantive_df))


# FIND WhETHER presence of TRAF2 help distinguish between endoctrine and structural

substantive_cells$TRAF2_status = ifelse(
  GetAssayData(substantive_cells, assay = "RNA", layer = "data")["TRAF2", ] > 0, 
  "TRAF2_positive", 
  "TRAF2_negative"
)

Idents(substantive_cells) = "TRAF2_status"

VlnPlot(substantive_cells, features = c("NRTN", "CAPRIN2"), pt.size = 0.1)





# VERIFYING IF NRTN NETWORK IS WHINE CHUCK SPECIFIC

confirmed_structural = c("NRTN", "SHISA8", "DDAH1", "WIF1", "ATP7B", "HIP1R", "PRKDC", "NUMB", "LOC140071013", "HMGCLL1") 

# Add a module score for network
substantive_cells = AddModuleScore(
  object = substantive_cells,
  features = list(confirmed_structural),
  name = "Verified_structural_score"
)

VlnPlot(substantive_cells, features = "Verified_structural_score1", group.by = "orig.ident", pt.size = 0) +
  labs(title = "Confirmed structural's network activation by sounds")



# Check dopamine now

drd2_expr = GetAssayData(wc_substantive, assay = "RNA", slot = "data")["DRD2", ]

# Categorize cells, if expressing (> 0), if not expressing (== 0)
wc_substantive$DRD2_status = ifelse(drd2_expr > 0, "DRD2_positive", "DRD2_negative")

Idents(wc_substantive) = "DRD2_status"

# Statistical test for growth and JNK
sensitivity_dop_stats = FindMarkers(
  wc_substantive,
  ident.1 = "DRD2_positive",
  ident.2 = "DRD2_negative",
  logfc.threshold = 0
)

dopamine_network = sensitivity_dop_stats %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
print(head(dopamine_network, 20))


VlnPlot(
  wc_substantive, 
  features = c("NRTN", "CAPRIN2"), 
  group.by = "DRD2_status", 
  pt.size = 0
)



# What fires NTRN?

wc_substantive$NRTN_Status = ifelse(
  GetAssayData(wc_substantive, assay = "RNA", layer = "data")["NRTN", ] > 0, 
  "NRTN_fired", 
  "NRTN_failed"
)

Idents(wc_substantive) = "NRTN_Status"
structural_triggers = FindMarkers(wc_substantive, ident.1 = "NRTN_fired", ident.2 = "NRTN_failed", logfc.threshold = 0.25)

nrtn_network = structural_triggers %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
print(head(nrtn_network, 20))



# Correlate nrtn network with all genes

structural_scores = substantive_cells$Verified_structural_score1
full_matrix = as.matrix(GetAssayData(substantive_cells, assay = "RNA", layer = "data"))

network_correlations = cor(structural_scores, t(full_matrix))

top_switches = sort(network_correlations[1, ], decreasing = TRUE)[2:31] 
print(as.data.frame(top_switches))


















# Fracture the yellow wgcna module

neural_modules = GetModules(neural_obj)

yellow_gene_names = neural_modules %>% 
  dplyr::filter(module == "yellow") %>%
  dplyr::pull(gene_name) # turns column into a text list


wc_substantive = subset(substantive_cells, subset = orig.ident == "Whine-chuck")
yellow_matrix = t(as.matrix(GetAssayData(wc_substantive, assay = "RNA", layer = "data")[yellow_gene_names, ]))

# Adjacency and Topological Overlap Matrix
# how tightly these genes share connections
adjacency_matrix = adjacency(yellow_matrix, type = "unsigned", power = 6) # Power 6 soft power
TOM = TOMsimilarity(adjacency_matrix)
dissTOM = 1 - TOM

# Hierarchical Clustering to group 138 genes into sub-communities
geneTree = hclust(as.dist(dissTOM), method = "average")

# Force tree to cut into sub-module
# minClusterSize forces it to find meaningful groups and not pairs of genes
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = 5)

# Building the table of exactly which genes fell into which fractured sub-module
yellow_fracture_table = data.frame(Gene = yellow_gene_names, Sub_community = dynamicMods)

# where our genes landed
print(yellow_fracture_table %>% filter(Gene %in% c("CAPRIN2", "NRTN")))


caprin2_network = yellow_fracture_table %>% filter(Sub_community == 0) %>% pull(Gene)
nrtn_network = yellow_fracture_table %>% filter(Sub_community == 1) %>% pull(Gene)

print(caprin2_network)
print(nrtn_network)





# Use those network to find receptors and makers


wc_substantive = AddModuleScore(wc_substantive, features = list(caprin2_network), name = "Endocrine_machine")
wc_substantive = AddModuleScore(wc_substantive, features = list(nrtn_network), name = "Structural_machine")

full_matrix = as.matrix(GetAssayData(wc_substantive, assay = "RNA", layer = "data"))

# ANTI CIRCULAR FILTER
upstream_matrix_endo = full_matrix[!(rownames(full_matrix) %in% caprin2_network), ]
upstream_matrix_struct = full_matrix[!(rownames(full_matrix) %in% nrtn_network), ]

endocrine_cor = cor(wc_substantive$Endocrine_machine1, t(upstream_matrix_endo))
structural_cor = cor(wc_substantive$Endocrine_machine1, t(upstream_matrix_struct))

top_endocrine_controllers = sort(endocrine_cor[1, ], decreasing = TRUE)[1:20]
top_structural_controllers = sort(structural_cor[1, ], decreasing = TRUE)[1:20]

endo_df = data.frame(tungara_symbol = names(top_endocrine_controllers), corr = top_endocrine_controllers)
struct_df = data.frame(tungara_symbol = names(top_structural_controllers), corr = top_structural_controllers)


print(head(endo_df, 15))
print(head(struct_df, 15))






# Confirming that the networks are only responding to Whine-chuck


substantive_cells = AddModuleScore(
  object = substantive_cells, 
  features = list(caprin2_network), 
  name = "True_endocrine_machine"
)

substantive_cells = AddModuleScore(
  object = substantive_cells, 
  features = list(nrtn_network), 
  name = "True_structural_machine"
)

VlnPlot(substantive_cells, 
        features = c("True_endocrine_machine1", "True_structural_machine1"), 
        group.by = "orig.ident", 
        pt.size = 0)

 
# Testing receptors that we found highly correlated and fold change in nrtn + some candidate genes


target_switches = c("GLP2R", "LTB4R2", "SP4", "SMAD2", "KDM6B")

VlnPlot(substantive_cells, 
        features = target_switches, 
        group.by = "orig.ident", 
        pt.size = .1, 
        ncol = 3)















# Running DE between Whine-chuck and Whine on substnative cells

Idents(substantive_cells) = "orig.ident"
cascade_markers = FindMarkers(substantive_cells, ident.1 = "Whine-chuck", ident.2 = "Whine", logfc.threshold = 0.25)

sig_cascade = cascade_markers %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(Gene = rownames(.))

# Categorize the genes by biological layer with nomenclature
sig_cascade = sig_cascade %>%
  mutate(Biological_Layer = case_when(
    grepl("R$|HTR|DRD|GRM|GAB|CHRN|ADGR|NMUR", Gene) ~ "Layer 1: Surface Receptors",
    grepl("^KCN|^CAC|^SCN|^SLC|^ATP", Gene) ~ "Layer 1b: Ion Channels / Transporters",
    grepl("K$|MAPK|TRAF|CAMK|PIK3|AKAP", Gene) ~ "Layer 2: Intracellular Kinases / Relays",
    grepl("TF|SMAD|SP|E2F|KDM|ZNF|HDAC", Gene) ~ "Layer 3: Nuclear Switches / Epigenetics",
    TRUE ~ "Layer 4: Downstream Effectors"
  ))

print(sig_cascade %>% filter(Biological_Layer %in% c("Layer 1: Surface Receptors", "Layer 1b: Ion Channels / Transporters")) %>% arrange(desc(avg_log2FC)) %>% select(Gene, avg_log2FC, pct.1, pct.2) %>% head(15))

print(sig_cascade %>% filter(Biological_Layer == "Layer 2: Intracellular Kinases / Relays") %>% arrange(desc(avg_log2FC)) %>% select(Gene, avg_log2FC, pct.1, pct.2) %>% head(15))

print(sig_cascade %>% filter(Biological_Layer == "Layer 3: Nuclear Switches / Epigenetics") %>% arrange(desc(avg_log2FC)) %>% select(Gene, avg_log2FC, pct.1, pct.2) %>% head(15))





# Isolate expression matrix for the active Whine-chuck cells
wc_substantive = subset(substantive_cells, subset = orig.ident == "Whine-chuck")
wc_matrix = GetAssayData(wc_substantive, assay = "RNA", layer = "data")

# cells fired main drivers
caprin2_positive_cells = colnames(wc_matrix)[wc_matrix["CAPRIN2", ] > 0]
nrtn_positive_cells = colnames(wc_matrix)[wc_matrix["NRTN", ] > 0]

# cooccurrence?
caprin2_co_occurrence = rowSums(wc_matrix[, caprin2_positive_cells] > 0) / length(caprin2_positive_cells)
nrtn_co_occurrence = rowSums(wc_matrix[, nrtn_positive_cells] > 0) / length(nrtn_positive_cells)

# baseline
baseline_matrix = GetAssayData(subset(substantive_cells, subset=orig.ident=="Whine"), assay="RNA", layer="data")
baseline_pct = rowSums(baseline_matrix > 0) / ncol(baseline_matrix)

# probability table
dependency_df = data.frame(
  Gene = rownames(wc_matrix),
  CAPRIN2_Dependency = caprin2_co_occurrence,
  NRTN_Dependency = nrtn_co_occurrence,
  Baseline_Noise = baseline_pct
)

# 60% of active cells required this gene)
true_triggers_caprin = dependency_df %>% filter(CAPRIN2_Dependency > 0.6 & Baseline_Noise < 0.2) %>% arrange(desc(CAPRIN2_Dependency))
true_triggers_nrtn = dependency_df %>% filter(NRTN_Dependency > 0.6 & Baseline_Noise < 0.2) %>% arrange(desc(NRTN_Dependency))

print(head(true_triggers_caprin, 15))

print(head(true_triggers_nrtn, 15))



# We need to know exactly which of those two neuron types is running the GLP2R -> GCK -> MAPK8 cascade, 
# and which one is sitting upstream of it.
# If they all have negative Log2FCs, the GABAergic neurons are the receiver
# If they have positive Log2FCs, Excitatory neurons are engine.

wc_substantive = subset(substantive_cells, subset = orig.ident == "Whine-chuck")

wc_substantive$Neuron_Class = ifelse(
  grepl("Excitatory", wc_substantive$refined_annotation), 
  "Excitatory_Class", 
  "GABAergic_Class"
)

# Expression of both
Idents(wc_substantive) = "Neuron_Class"
class_markers = FindMarkers(wc_substantive, ident.1 = "Excitatory_Class", ident.2 = "GABAergic_Class", logfc.threshold = 0.25)

# Filter for significance and categorize
sig_class = class_markers %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(Gene = rownames(.)) %>%
  arrange(desc(avg_log2FC))

print(head(sig_class, 15))

print(tail(sig_class, 15))


# 1. Search for ligands
ligand_genes = grep("^GCG|^GLP|^SCT", rownames(neural_obj), value = TRUE)
print(ligand_genes)

if (length(ligand_genes) > 0) {
  # expression just for Whine-chuck 
  wc_global = subset(neural_obj, subset = orig.ident == "Whine-chuck")
  
  # average expression of the ligands across all cell clusters
  ligand_production = AverageExpression(wc_global, features = ligand_genes, group.by = "refined_annotation")$RNA
  
  ligand_df = as.data.frame(ligand_production)
  print(ligand_df)
} else {
  print("NO LIGANDS FOUND IN MATRIX ")
}




# look in neural_obj for GCG and SCT

wc_global = subset(neural_obj, subset = orig.ident == "Whine-chuck")

# ligands for receptors
target_ligands = c("GCG", "SCT") 
valid_ligands = target_ligands[target_ligands %in% rownames(wc_global)]

# average expression
ligand_expr = AverageExpression(wc_global, features = valid_ligands, group.by = "refined_annotation")$RNA

ligand_df = as.data.frame(ligand_expr)
ligand_df$Ligand = rownames(ligand_df)
ligand_long = pivot_longer(ligand_df, cols = -Ligand, names_to = "Cell_Type", values_to = "Expression")

ligand_long$Cell_Type = gsub("\\.", " ", ligand_long$Cell_Type)

top_senders = ligand_long %>%
  filter(Expression > 0) %>%
  arrange(desc(Expression))

print(head(top_senders, 15))



# #1 ligand producing cell type
top_sender_name = top_senders$Cell_Type[1]
print(paste("SENDER:", top_sender_name))

# Sender cluster across all sound treatments
sender_cells = subset(neural_obj, subset = refined_annotation == top_sender_name)
Idents(sender_cells) = "orig.ident"

# De on Sender
sender_response = FindMarkers(sender_cells, ident.1 = "Whine-chuck", ident.2 = "Whine", logfc.threshold = 0.25)

# Filter for significance
sig_sender = sender_response %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(Gene = rownames(.)) %>%
  arrange(desc(avg_log2FC))

print(head(sig_sender, 20))





# Whole-Brain Ligand findin

global_wc = subset(tunga_frog_merged_blast, subset = orig.ident == "Whine-chuck")

# target ligands
target_ligands = c("GCG", "SCT") 
valid_ligands = target_ligands[target_ligands %in% rownames(global_wc)]

global_ligand_expr = AverageExpression(global_wc, features = valid_ligands, group.by = "refined_annotation")$RNA

global_ligand_df = as.data.frame(global_ligand_expr)
global_ligand_df$Ligand = rownames(global_ligand_df)
global_ligand_long = pivot_longer(global_ligand_df, cols = -Ligand, names_to = "Cell_Type", values_to = "Expression")

global_ligand_long$Cell_Type = gsub("\\.", " ", global_ligand_long$Cell_Type)

top_global_senders = global_ligand_long %>%
  filter(Expression > 0) %>%
  arrange(desc(Expression))

print(head(top_global_senders, 15))


wc_substantive = subset(substantive_cells, subset = orig.ident == "Whine-chuck")

# Is AVP exclusive to Whine-chuck response
# TAC1 and SST to see what are they doing
VlnPlot(substantive_cells, 
        features = c("CAPRIN2", "AVP", "TAC1", "SST"), 
        group.by = "orig.ident", 
        pt.size = 0, ncol = 4)

# Does AVP require the Endocrine?
wc_matrix = GetAssayData(wc_substantive, assay = "RNA", layer = "data")

# Find cells successfully building the AVP hormone
avp_positive_cells = colnames(wc_matrix)[wc_matrix["AVP", ] > 0]

# Calculate what percentage of those AVP+ cells are running CAPRIN2
if(length(avp_positive_cells) > 0) {
  caprin2_avp_overlap = sum(wc_matrix["CAPRIN2", avp_positive_cells] > 0) / length(avp_positive_cells)
  print(paste("AVP CELLS BY CAPRIN2 --- :", round(caprin2_avp_overlap * 100, 2), "%"))
} else {
  print("AVP NOT DETECTED")
}
at this entire cascade exists solely to pump Vasopressin.


Idents(substantive_cells) = "orig.ident"
cascade_markers = FindMarkers(substantive_cells, ident.1 = "Whine-chuck", ident.2 = "Whine", logfc.threshold = 0.25)

wnt_blueprint = cascade_markers %>%
  mutate(Gene = rownames(.)) %>%
  filter(grepl("^WNT|^FZD|^LRP|^DVL|^CTNNB|^WIF|^DKK", Gene)) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  select(Gene, avg_log2FC, pct.1, pct.2)

print(wnt_blueprint)

wc_matrix = GetAssayData(wc_substantive, assay = "RNA", layer = "data")

proposed_cascade = c("GLP2R", "GCK", "MAPK8", "KDM6B", "CAPRIN2", "NRTN")

valid_cascade = proposed_cascade[proposed_cascade %in% rownames(wc_matrix)]

dense_subset = as.matrix(wc_matrix[valid_cascade, , drop = FALSE])

binary_matrix = (dense_subset > 0) * 1 

prob_matrix = matrix(0, nrow = length(valid_cascade), ncol = length(valid_cascade))
rownames(prob_matrix) = valid_cascade
colnames(prob_matrix) = valid_cascade

for (row_gene in valid_cascade) {
  for (col_gene in valid_cascade) {
    row_active_cells = which(binary_matrix[row_gene, ] == 1)
    
    if (length(row_active_cells) > 0) {
      co_active = sum(binary_matrix[col_gene, row_active_cells] == 1)
      prob_matrix[row_gene, col_gene] = co_active / length(row_active_cells)
    }
  }
}

print(round(prob_matrix, 3))

pheatmap(prob_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = TRUE,
         main = "Cascade genes checking")

caprin2_network_genes = yellow_fracture_table %>% 
  filter(Sub_community == 0) %>% 
  pull(Gene)

caprin2_expr = as.numeric(wc_matrix["CAPRIN2", ])
network_matrix_c = as.matrix(wc_matrix[rownames(wc_matrix) %in% caprin2_network_genes, , drop=FALSE])

caprin_correlations = cor(caprin2_expr, t(network_matrix_c))

caprin_downstream = data.frame(
  Gene = colnames(caprin_correlations),
  CAPRIN2_Corr = as.vector(caprin_correlations)
) %>% arrange(desc(CAPRIN2_Corr))

caprin_downstream = caprin_downstream %>%
  mutate(Payload_Type = case_when(
    grepl("^RPL|^RPS|^EIF", Gene) ~ "Ribosomal/Translation Engine",
    grepl("POMC|OXT|AVP|CRH|TRH|NPY|VIP|TAC", Gene) ~ "Neuropeptide/Hormone",
    TRUE ~ "Other"
  ))

print(head(caprin_downstream %>% filter(Payload_Type != "Other"), 25))

nrtn_network_genes = yellow_fracture_table %>% 
  filter(Sub_community == 1) %>% 
  pull(Gene)

nrtn_expr = as.numeric(wc_matrix["NRTN", ])
network_matrix_n = as.matrix(wc_matrix[rownames(wc_matrix) %in% nrtn_network_genes, , drop=FALSE])

nrtn_correlations = cor(nrtn_expr, t(network_matrix_n))

nrtn_downstream = data.frame(
  Gene = colnames(nrtn_correlations),
  NRTN_Corr = as.vector(nrtn_correlations)
) %>% arrange(desc(NRTN_Corr))

nrtn_downstream = nrtn_downstream %>%
  mutate(Structure_Type = case_when(
    grepl("^ACT|^TUB|^MYO", Gene) ~ "Cytoskeleton (Actin/Tubulin)",
    grepl("^SYP|^SST|^VAMP|^SNAP", Gene) ~ "Synaptic Vesicle/Fusion",
    grepl("^NLGN|^NRXN|^DLG", Gene) ~ "Synapse Scaffolding",
    TRUE ~ "Other"
  ))

print(head(nrtn_downstream %>% filter(Structure_Type != "Other"), 25))

caprin_downstream_unfiltered = data.frame(
  tungara_symbol = colnames(caprin_correlations),
  CAPRIN2_Corr = as.vector(caprin_correlations)
) %>% arrange(desc(CAPRIN2_Corr))

print(head(caprin_downstream_unfiltered, 25))

translated_results10 = left_join(caprin_downstream_unfiltered, translation_table, by = "tungara_symbol")

final_table = translated_results10 %>%
  select(tungara_symbol, xenopus_symbol, CAPRIN2_Corr)

print(head(final_table, 30))

nrtn_downstream_unfiltered = data.frame(
  tungara_symbol = colnames(nrtn_correlations),
  NRTN_Corr = as.vector(nrtn_correlations)
) %>% arrange(desc(NRTN_Corr))

print(head(nrtn_downstream_unfiltered, 25))

translated_results11 = left_join(nrtn_downstream_unfiltered, translation_table, by = "tungara_symbol")

final_table = translated_results11 %>%
  select(tungara_symbol, xenopus_symbol, NRTN_Corr)

print(head(final_table, 30))

wc_matrix = as.matrix(GetAssayData(wc_substantive, assay = "RNA", layer = "data"))
caprin2_expr = as.numeric(wc_matrix["CAPRIN2", ])

all_caprin2_cor = cor(caprin2_expr, t(wc_matrix))

caprin2_true_partners = data.frame(
  Gene = colnames(all_caprin2_cor),
  Correlation = as.vector(all_caprin2_cor)
) %>% filter(Gene != "CAPRIN2") %>% arrange(desc(Correlation))

print(head(caprin2_true_partners, 25))

print(tail(caprin2_true_partners, 15))

nrtn_expr = as.numeric(wc_matrix["NRTN", ])

co_fire_table = data.frame(
  Cell = colnames(wc_matrix),
  CAPRIN2_Fired = caprin2_expr > 0,
  NRTN_Fired = nrtn_expr > 0,
  Cell_Type = wc_substantive$refined_annotation 
)

print(table(CAPRIN2_Active = co_fire_table$CAPRIN2_Fired, NRTN_Active = co_fire_table$NRTN_Fired))

print(table(co_fire_table$Cell_Type[co_fire_table$NRTN_Fired]))

print(table(co_fire_table$Cell_Type[co_fire_table$CAPRIN2_Fired]))

all_nrtn_cor = cor(nrtn_expr, t(wc_matrix))

switch_df = data.frame(
  Gene = colnames(all_nrtn_cor),
  Correlation = as.vector(all_nrtn_cor)
) %>% arrange(desc(Correlation))

true_switches = switch_df %>%
  filter(grepl("TF$|ZNF|KDM|HDAC|STAT|SMAD|MAPK|CAMK|CREB|FOS|JUN|EGR", Gene))

print(head(true_switches, 20))

true_receptors = switch_df %>%
  filter(grepl("R$|HTR|DRD|GRM|GAB|CHRN|ADGR", Gene)) %>%
  filter(!grepl("RPL|RPS|TAR", Gene)) 

print(head(true_receptors, 20))

gaba_engine = subset(wc_substantive, subset = refined_annotation == "GABAergic neuron_tac1 high")
gaba_matrix = as.matrix(GetAssayData(gaba_engine, assay = "RNA", layer = "data"))

real_cascade = c("GLP2R", "HIP1R", "CAMK1", "MAPK8", "CREB1", "ZNF292", "NRTN", "DDAH1", "SLC22A2")
valid_real_cascade = real_cascade[real_cascade %in% rownames(gaba_matrix)]

cascade_cor_matrix = cor(t(gaba_matrix[valid_real_cascade, ]))

print(round(cascade_cor_matrix, 3))

nrtn_status_gaba = ifelse(gaba_matrix["NRTN", ] > 0, "NRTN_Active", "NRTN_Inactive")
gaba_engine = AddMetaData(gaba_engine, metadata = nrtn_status_gaba, col.name = "NRTN_State")

Idents(gaba_engine) = "NRTN_State"
nrtn_effect = FindMarkers(gaba_engine, ident.1 = "NRTN_Active", ident.2 = "NRTN_Inactive", features = c("SLC22A2", "DDAH1", "HIP1R", "GLP2R", "MAPK8"), logfc.threshold = 0)

print(nrtn_effect)

wc_global = subset(neural_obj, subset = orig.ident == "Whine-chuck")

dopamine_machinery = c("TH", "DDC", "DBH")
valid_da_genes = dopamine_machinery[dopamine_machinery %in% rownames(wc_global)]

da_expr = AverageExpression(wc_global, features = valid_da_genes, group.by = "refined_annotation")$RNA

da_df = as.data.frame(da_expr)
da_df$Gene = rownames(da_df)
da_long = pivot_longer(da_df, cols = -Gene, names_to = "Cell_Type", values_to = "Expression")

da_long$Cell_Type = gsub("\\.", " ", da_long$Cell_Type)

top_da_senders = da_long %>%
  group_by(Cell_Type) %>%
  summarize(Total_DA_Machinery = sum(Expression)) %>%
  arrange(desc(Total_DA_Machinery)) %>%
  filter(Total_DA_Machinery > 0)

print(head(top_da_senders, 15))

hormone_payloads = c("AVP", "OXT", "GNRH1", "CRH", "TRH", "POMC")
valid_hormones = hormone_payloads[hormone_payloads %in% rownames(wc_global)]

hormone_expr = AverageExpression(wc_global, features = valid_hormones, group.by = "refined_annotation")$RNA

hormone_df = as.data.frame(hormone_expr)
hormone_df$Gene = rownames(hormone_df)
hormone_long = pivot_longer(hormone_df, cols = -Gene, names_to = "Cell_Type", values_to = "Expression")

hormone_long$Cell_Type = gsub("\\.", " ", hormone_long$Cell_Type)

top_hormone_factories = hormone_long %>%
  group_by(Cell_Type) %>%
  summarize(Total_Hormone_Output = sum(Expression)) %>%
  arrange(desc(Total_Hormone_Output)) %>%
  filter(Total_Hormone_Output > 0)

print(head(top_hormone_factories, 15))

receptors_to_check = grep("^GABR|^GRIN|^GRIA|^HTR|^DRD", rownames(gaba_engine), value = TRUE)

receptor_shift = FindMarkers(gaba_engine, 
                              ident.1 = "NRTN_Active", 
                              ident.2 = "NRTN_Inactive", 
                              features = receptors_to_check, 
                              logfc.threshold = 0.1)

swallowed_receptors = receptor_shift %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC)

print(head(swallowed_receptors, 15))

nrtn_down = FindMarkers(gaba_engine, 
                         ident.1 = "NRTN_Active", 
                         ident.2 = "NRTN_Inactive", 
                         logfc.threshold = 0.25)

swallowed_unbiased = nrtn_down %>%
  filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
  arrange(avg_log2FC) %>%
  mutate(tungara_symbol = rownames(.))

print(head(swallowed_unbiased %>% select(tungara_symbol, avg_log2FC, pct.1, pct.2), 20))

translated_results12 = left_join(swallowed_unbiased, translation_table, by = "tungara_symbol")

final_table = translated_results12 %>%
  select(tungara_symbol, xenopus_symbol, avg_log2FC, pct.1, pct.2)

print(head(final_table, 30))

micro_circuit = subset(wc_substantive, subset = refined_annotation %in% c("GABAergic neuron_sst high", "GABAergic neuron_tac1 high"))
Idents(micro_circuit) = "refined_annotation"

handshake = FindMarkers(micro_circuit, 
                         ident.1 = "GABAergic neuron_sst high", 
                         ident.2 = "GABAergic neuron_tac1 high", 
                         logfc.threshold = 0.25)

handshake = handshake %>% mutate(Gene = rownames(.))

sst_makers = handshake %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  arrange(desc(avg_log2FC))

tac1_makers = handshake %>% 
  filter(p_val_adj < 0.05 & avg_log2FC < 0) %>% 
  arrange(avg_log2FC)

print(head(sst_makers %>% select(Gene, avg_log2FC, pct.1, pct.2), 15))

print(head(tac1_makers %>% select(Gene, avg_log2FC, pct.1, pct.2), 15))

`%notin%` = Negate(`%in%`)

Idents(wc_substantive) = "Epicenter"
rest_of_brain = subset(wc_global, subset = refined_annotation %notin% c("GABAergic neuron_sst high", "GABAergic neuron_tac1 high", "Excitatory neuron_tac1 high"))
Idents(rest_of_brain) = "Other"

combined_temp = merge(wc_substantive, y = rest_of_brain)

combined_temp = JoinLayers(combined_temp)

master_markers = FindMarkers(combined_temp, ident.1 = "Epicenter", ident.2 = "Other", logfc.threshold = 0.5)

master_switches = master_markers %>%
  mutate(Gene = rownames(.)) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  filter(grepl("TF$|ZNF|KDM|HDAC|STAT|SMAD|MAPK|CAMK|CREB|FOS|JUN|EGR|PAX|GATA", Gene)) %>%
  arrange(desc(avg_log2FC))

print(head(master_switches %>% select(Gene, avg_log2FC, pct.1, pct.2), 20))

fixed_plasticity_targets = c("SLC22A2", "LOC140108393", "DDAH1", "HIP1R")

fig_5_fixed = VlnPlot(gaba_engine, 
                       features = fixed_plasticity_targets, 
                       group.by = "NRTN_State", 
                       pt.size = 0, 
                       ncol = 4) +
  labs(title = "Fixed: chemical plasticity",
       subtitle = "LOC140108393 is ITGAD. Watch it drop.")

print(fig_5_fixed)

nrtn_triggers = FindMarkers(gaba_engine, 
                             ident.1 = "NRTN_Fired", 
                             ident.2 = "NRTN_Failed", 
                             logfc.threshold = 0.25,
                             min.pct = 0.25) 

dense_receptors = nrtn_triggers %>%
  mutate(Gene = rownames(.)) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  filter(grepl("R$|HTR|DRD|GRM|GAB|CHRN|ADGR", Gene)) %>%
  filter(!grepl("RPL|RPS", Gene)) %>% 
  arrange(desc(pct.1)) 

print(head(dense_receptors %>% select(Gene, avg_log2FC, pct.1, pct.2), 10))

sst_cells = subset(wc_substantive, subset = refined_annotation == "GABAergic neuron_sst high")
Idents(sst_cells) = "orig.ident"

sst_cells$

sst_payload = FindMarkers(sst_cells, 
                           ident.1 = "Whine-chuck", 
                           ident.2 = "Whine", 
                           logfc.threshold = 0.25,
                           min.pct = 0.3) 

true_sst_makers = sst_payload %>%
  mutate(Gene = rownames(.)) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC))

print(head(true_sst_makers %>% select(Gene, avg_log2FC, pct.1, pct.2), 15))

Idents(substantive_cells) = "refined_annotation"
sst_all_treatments = subset(substantive_cells, idents = "GABAergic neuron_sst high")

Idents(sst_all_treatments) = "orig.ident"

sst_payload = FindMarkers(sst_all_treatments, 
                           ident.1 = "Whine-chuck", 
                           ident.2 = "Whine", 
                           logfc.threshold = 0.25,
                           min.pct = 0.25) 

true_sst_makers = sst_payload %>%
  mutate(Gene = rownames(.)) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC))

print(head(true_sst_makers %>% select(Gene, avg_log2FC, pct.1, pct.2), 15))

memory_wire_genes_fixed = c("ADGRA3", "CAMK1", "CREB1", "MAPK8", "ZNF292")

fig_4_fixed = DotPlot(gaba_engine, features = memory_wire_genes_fixed, group.by = "orig.ident") +
  coord_flip() +
  labs(title = "internal wire",
       subtitle = "ADGRA3 triggers Calcium/CaMK-CREB-MAPK")

print(fig_4_fixed)

plasticity_targets_fixed = c("SLC22A2", "DDAH1", "HIP1R")

fig_5_fixed = VlnPlot(gaba_engine, 
                       features = plasticity_targets_fixed, 
                       group.by = "NRTN_State", 
                       pt.size = 0, 
                       ncol = 3) +
  labs(title = "Chemical plasticity lockdown")

print(fig_5_fixed)

tac1_all_treatments = subset(substantive_cells, idents = "GABAergic neuron_tac1 high")
Idents(tac1_all_treatments) = "orig.ident"

tac1_response = FindMarkers(tac1_all_treatments, 
                             ident.1 = "Whine-chuck", 
                             ident.2 = "Whine", 
                             logfc.threshold = 0.25)

significant_tac1 = tac1_response %>%
  mutate(Gene = rownames(.)) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC))

print(head(significant_tac1 %>% select(Gene, avg_log2FC, pct.1, pct.2), 15))

micro_circuit_all = subset(substantive_cells, subset = refined_annotation %in% c("GABAergic neuron_sst high", "GABAergic neuron_tac1 high"))

micro_circuit_all$Condition_Cell = paste(micro_circuit_all$refined_annotation, micro_circuit_all$orig.ident, sep = "_")
Idents(micro_circuit_all) = "Condition_Cell"

true_engine_genes = c("CAPRIN2", "NRTN", "PCSK6", "SHISA8", "HIP1R")

fig_sst_explosion = DotPlot(micro_circuit_all, features = true_engine_genes) + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "CAPRIN2 and NRTN only explode in SST cells during Whine-chuck")

print(fig_sst_explosion)

excitatory_cells = subset(substantive_cells, idents = "Excitatory neuron_tac1 high")
Idents(excitatory_cells) = "orig.ident"

excitatory_response = FindMarkers(excitatory_cells, 
                                   ident.1 = "Whine-chuck", 
                                   ident.2 = "Whine", 
                                   logfc.threshold = 0.25)

sig_excitatory = excitatory_response %>%
  mutate(Gene = rownames(.)) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC))

print(head(sig_excitatory %>% select(Gene, avg_log2FC, pct.1, pct.2), 20))

shared_engine_genes = c("CAPRIN2", "NRTN", "HIP1R", "WIF1", "SLC22A2")

print(sig_excitatory %>% 
        filter(Gene %in% shared_engine_genes) %>% 
        select(Gene, avg_log2FC, pct.1, pct.2))

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

Idents(substantive_cells) = "orig.ident"

wc_epicenter = subset(substantive_cells, subset = orig.ident == "Whine-chuck")
wc_epicenter$CellClass = ifelse(grepl("Excitatory", wc_epicenter$refined_annotation), "Excitatory", "Inhibitory")
Idents(wc_epicenter) = "CellClass"

exit_wire_markers = FindMarkers(wc_epicenter, ident.1 = "Excitatory", ident.2 = "Inhibitory", logfc.threshold = 0.25)

projection_genes = exit_wire_markers %>%
  mutate(Gene = rownames(.)) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  filter(grepl("^EPH|^SEMA|^SLIT|^ROBO|^PLXN|^NTRK|^NRXN|^NLGN|^SYT|^VAMP", Gene)) %>%
  arrange(desc(avg_log2FC))

print(head(projection_genes %>% select(Gene, avg_log2FC, pct.1, pct.2), 20))

Idents(substantive_cells) = "orig.ident"
whine_baseline = FindMarkers(substantive_cells, 
                              ident.1 = "Whine", 
                              ident.2 = "Whine-chuck", 
                              logfc.threshold = 0.25)

sig_whine = whine_baseline %>%
  mutate(tungara_symbol = rownames(.)) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC))

print(head(sig_whine %>% select(tungara_symbol, avg_log2FC, pct.1, pct.2), 20))

translated_results14 = left_join(head(sig_whine %>% select(tungara_symbol, avg_log2FC, pct.1, pct.2)), translation_table, by = "tungara_symbol")

final_table = translated_results14 %>% select(tungara_symbol, xenopus_symbol, avg_log2FC, pct.1, pct.2) %>%
  arrange(desc(avg_log2FC))

print(head(final_table, 30))

known_ligands = c("FLRT1", "FLRT2", "FLRT3", "GCG", "SCT", "BDNF", "NTF3", "NTF4", "NGF", "WNT5A", "WNT7A")
valid_ligands = known_ligands[known_ligands %in% rownames(wc_global)]

if(length(valid_ligands) > 0) {
  ligand_expr = AverageExpression(wc_global, features = valid_ligands, group.by = "refined_annotation")$RNA
  ligand_df = as.data.frame(ligand_expr)
  ligand_df$Ligand = rownames(ligand_df)
  ligand_long = pivot_longer(ligand_df, cols = -Ligand, names_to = "Cell_Type", values_to = "Expression") %>% filter(Expression > 0)
  
  print(ligand_long %>% arrange(desc(Expression)) %>% head(20))
}

epicenter_matrix = as.matrix(GetAssayData(wc_epicenter, assay = "RNA", layer = "data"))
caprin_epicenter_expr = as.numeric(epicenter_matrix["CAPRIN2", ])

caprin_cor = cor(caprin_epicenter_expr, t(epicenter_matrix))
caprin_secretome = data.frame(Gene = colnames(caprin_cor), Corr = as.vector(caprin_cor)) %>%
  filter(grepl("^FGF|^VEGF|^IGF|^TGF|^BMP|^WNT|^CXCL|^CCL|^IL|^NPY|^PENK|^PDYN", Gene)) %>%
  arrange(desc(Corr))

print(head(caprin_secretome, 20))

nrtn_epicenter_expr = as.numeric(epicenter_matrix["NRTN", ])

engine_score = caprin_epicenter_expr + nrtn_epicenter_expr
engine_cor = cor(engine_score, t(epicenter_matrix))

master_tf = data.frame(Gene = colnames(engine_cor), Corr = as.vector(engine_cor)) %>%
  filter(grepl("TF$|ZNF|KDM|HDAC|STAT|SMAD|MAPK|CAMK|CREB|FOS|JUN|EGR|PAX|GATA", Gene)) %>%
  arrange(desc(Corr))

print(head(master_tf, 20))

library(patchwork)

Idents(wc_global) = "refined_annotation"

spatial_glyph_genes = c("PAX2", "PAX3", "PAX5", "PAX7", "GATA2", "GATA3")

fig_1 = DotPlot(wc_global, features = spatial_glyph_genes) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "epicenter in the Torus Semicircularis")

print(fig_1)

fig_2 = VlnPlot(wc_global, 
                 features = "SLC6A5", 
                 pt.size = 1, 
                 sort = TRUE) +
  theme(legend.position = "none") +
  labs(title = "SLC6A5 expression, Whine-chuck from green noise")

print(fig_2)

micro_circuit = subset(wc_global, subset = refined_annotation %in% c("GABAergic neuron_sst high", "GABAergic neuron_tac1 high"))
Idents(micro_circuit) = "refined_annotation"

factory_genes = c("EGR3", "EGR4", "TH", "DDC", "POMC", "PLAT")

fig_3 = DotPlot(micro_circuit, features = factory_genes) + 
  coord_flip() +
  labs(title = "SST cells")

print(fig_3)

gaba_engine = subset(wc_global, subset = refined_annotation == "GABAergic neuron_tac1 high")

memory_wire_genes = c("GLP2R", "CAMK1", "CREB1", "MAPK8", "ZNF292")

fig_4 = DotPlot(gaba_engine, features = memory_wire_genes, group.by = "orig.ident") +
  coord_flip() +
  labs(title = "CaMK-CREB-MAPK memory")

print(fig_4)

gaba_engine$NRTN_State = ifelse(
  GetAssayData(gaba_engine, assay = "RNA", layer = "data")["NRTN", ] > 0, 
  "NRTN_Fired", 
  "NRTN_Failed"
)
Idents(gaba_engine) = "NRTN_State"
plasticity_targets = c("SLC22A2", "LOC140108393", "DDAH1", "HIP1R")

fig_5 = VlnPlot(gaba_engine, 
                 features = plasticity_targets, 
                 group.by = "NRTN_State", 
                 pt.size = 0, 
                 ncol = 2) +
  labs(title = "NRTN traps dopamine")

print(fig_5)

epicenter_all = subset(substantive_cells, subset = refined_annotation %in% c("GABAergic neuron_sst high", "GABAergic neuron_tac1 high", "Excitatory neuron_tac1 high"))

epicenter_wc = subset(epicenter_all, subset = orig.ident == "Whine-chuck")
dense_epi_wc = as.matrix(GetAssayData(epicenter_wc, assay = "RNA", layer = "data"))

wire_genes = c("ADGRA3", "CAMK1", "CREB1", "MAPK8", "ZNF292", "NRTN")
valid_wire = wire_genes[wire_genes %in% rownames(dense_epi_wc)]

wire_cor = cor(t(dense_epi_wc[valid_wire, ]))

print(round(wire_cor, 3))

nrtn_status = ifelse(dense_epi_wc["NRTN", ] > 0, "NRTN_Active", "NRTN_Inactive")
epicenter_wc = AddMetaData(epicenter_wc, metadata = nrtn_status, col.name = "Global_NRTN_State")

Idents(epicenter_wc) = "Global_NRTN_State"

lockdown_targets = c("SLC22A2", "HIP1R", "DDAH1", "SYT7", "EPHA8", "WIF1")
valid_targets = lockdown_targets[lockdown_targets %in% rownames(epicenter_wc)]

universal_phenotype = FindMarkers(epicenter_wc, 
                                   ident.1 = "NRTN_Active", 
                                   ident.2 = "NRTN_Inactive", 
                                   features = valid_targets, 
                                   logfc.threshold = 0)

print(universal_phenotype %>% select(avg_log2FC, pct.1, pct.2, p_val_adj))

Idents(epicenter_all) = "orig.ident"

baseline_test = FindMarkers(epicenter_all, 
                             ident.1 = "Whine", 
                             ident.2 = "Whine-chuck", 
                             features = c("BOC", "LOC140111485"), 
                             logfc.threshold = 0)

print(baseline_test %>% select(avg_log2FC, pct.1, pct.2, p_val_adj))

whole_brain_all = tunga_frog_merged_blast

whole_brain_all$is_complex_sound = whole_brain_all$orig.ident != "Whine"
whole_brain_all$sound_state = ifelse(whole_brain_all$is_complex_sound, "Complex", "Baseline")
Idents(whole_brain_all) = "sound_state"

complex_brain = subset(whole_brain_all, idents = "Complex")

dopamine_genes = c("TH", "DDC", "DBH")
dopamine_expression = AverageExpression(complex_brain, features = dopamine_genes, group.by = "refined_annotation")$RNA
dopamine_table = as.data.frame(dopamine_expression)
dopamine_table$Gene = rownames(dopamine_table)
print(head(dopamine_table, 20))

ligand_genes = c("GCG", "SCT", "FLRT1", "BDNF", "NTF3")
ligand_expression = AverageExpression(complex_brain, features = ligand_genes, group.by = "refined_annotation")$RNA
ligand_table = as.data.frame(ligand_expression)
ligand_table$Gene = rownames(ligand_table)
print(head(ligand_table, 20))

receiver_cells = subset(complex_brain, subset = refined_annotation == "GABAergic neuron_tac1 high")
receiver_matrix = as.matrix(GetAssayData(receiver_cells, assay = "RNA", layer = "data"))
caprin2_vector = as.numeric(receiver_matrix["CAPRIN2", ])

transport_genes = c("KIF5A", "DYNC1H1", "FMR1", "PURA", "ELAVL1", "ELAVL2")
valid_transport = transport_genes[transport_genes %in% rownames(receiver_matrix)]
transport_cor = cor(caprin2_vector, t(receiver_matrix[valid_transport, , drop = FALSE]))
transport_result = data.frame(Gene = colnames(transport_cor), Score = as.vector(transport_cor))
print(transport_result)

granule_genes = c("G3BP1", "TIA1", "PABPC1", "CAPRIN1", "ATXN2", "YBX1")
valid_granules = granule_genes[granule_genes %in% rownames(receiver_matrix)]
granule_cor = cor(caprin2_vector, t(receiver_matrix[valid_granules, , drop = FALSE]))
granule_result = data.frame(Gene = colnames(granule_cor), Score = as.vector(granule_cor))
print(granule_result)

wnt_genes = c("CTNNB1", "GSK3B", "LRP5", "LRP6", "DVL1", "AXIN1")
valid_wnt = wnt_genes[wnt_genes %in% rownames(receiver_matrix)]
wnt_cor = cor(caprin2_vector, t(receiver_matrix[valid_wnt, , drop = FALSE]))
wnt_result = data.frame(Gene = colnames(wnt_cor), Score = as.vector(wnt_cor))
print(wnt_result)

glia_signal_genes = c("CX3CL1", "CD200", "SIRPA", "CSF1", "IL34")
valid_glia = glia_signal_genes[glia_signal_genes %in% rownames(receiver_matrix)]
glia_cor = cor(caprin2_vector, t(receiver_matrix[valid_glia, , drop = FALSE]))
glia_result = data.frame(Gene = colnames(glia_cor), Score = as.vector(glia_cor))
print(glia_result)

receiver_matrix = as.matrix(GetAssayData(receiver_cells, assay = "RNA", layer = "data"))
caprin2_vector = as.numeric(receiver_matrix["CAPRIN2", ])

global_correlations = apply(receiver_matrix, 1, function(x) {
  if(sd(x) == 0) return(NA)
  return(cor(caprin2_vector, x))
})

global_correlations = na.omit(global_correlations)
sorted_correlations = sort(global_correlations, decreasing = TRUE)

top_30_partners = sorted_correlations[2:31]
dragnet_results = data.frame(Gene = names(top_30_partners), Correlation_Score = as.numeric(top_30_partners))

print(dragnet_results)

receiver_cells$caprin2_level = caprin2_vector

caprin2_threshold = quantile(caprin2_vector[caprin2_vector > 0], 0.75)

receiver_cells$activation_state = ifelse(receiver_cells$caprin2_level >= caprin2_threshold, "High_Saturation",
                                         ifelse(receiver_cells$caprin2_level > 0, "Low_Saturation", "Baseline_Off"))

Idents(receiver_cells) = "activation_state"

caprin2_cascade = FindMarkers(receiver_cells, 
                              ident.1 = "High_Saturation", 
                              ident.2 = "Low_Saturation", 
                              logfc.threshold = 0.5, 
                              only.pos = TRUE)

print(head(caprin2_cascade, 20))