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




# CHECKPOINT 5.5: Annotation
# ==============================================================================

if (file.exists("tunga_frog_annotated.rds")) {
  
  tunga_frog_merged = readRDS("tunga_frog_annotated.rds")
  
} else {
  
  tunga_frog_merged = readRDS("tunga_frog_clustered.rds")
  
  brain = read_excel("xenopus_brain_markers.xlsx", sheet = "Brain")
  brain = brain %>% fill(`Cell type`) 
  
  orthos = read.csv("xenopus_gene_to_tungara_gene.csv")
  brain_mapped = brain %>%
    left_join(
      orthos,
      by = c("gene" = "xenopus_gene"),
      relationship = "many-to-many"
    ) %>%
    filter(!is.na(tungara_gene)) %>%
    distinct(`Cell type`, gene, tungara_gene, .keep_all = TRUE)
  
  pituitary_keywords = "Gonadotroph|Thyrotroph|Endocrine|Growth hormone|Prolactin|Melanotrope"
  
  brain_mapped_strong = brain_mapped %>%
    filter(!grepl(pituitary_keywords, `Cell type`, ignore.case = TRUE)) %>%
    group_by(`Cell type`) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 50)
  
  sctype_db = brain_mapped_strong %>%
    group_by(`Cell type`) %>%
    summarize(
      geneSymbolmore1 = paste(unique(tungara_gene), collapse = ","),
      .groups = "drop"
    ) %>%
    mutate(
      tissueType = "Brain",
      geneSymbolmore2 = ""
    ) %>%
    rename(cellName = `Cell type`) %>%
    select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2)
  
  openxlsx::write.xlsx(sctype_db, "tungara_brain_sctype_markers_clean.xlsx", rowNames = FALSE)
  
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  gs_list = gene_sets_prepare(
    "tungara_brain_sctype_markers_clean.xlsx",
    "Brain"
  )
  
  es.max = sctype_score(
    scRNAseqData = tunga_frog_merged[["RNA"]]$scale.data,
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  
  sctype_top = do.call("rbind", lapply(unique(tunga_frog_merged$seurat_clusters), function(cl) {
    cells = rownames(tunga_frog_merged@meta.data)[tunga_frog_merged$seurat_clusters == cl]
    scores = sort(rowSums(es.max[, cells, drop = FALSE]), decreasing = TRUE)
    
    data.frame(
      cluster = cl,
      cell_type = names(scores)[1],
      score = scores[1]
    )
  }))
  
  tunga_frog_merged$sctype_label = sctype_top$cell_type[
    match(tunga_frog_merged$seurat_clusters, sctype_top$cluster)
  ]
  
  tunga_frog_merged@meta.data = tunga_frog_merged@meta.data %>%
    mutate(broad_annotation = case_when(
      grepl("Excitatory neuron", sctype_label, ignore.case = TRUE) ~ "Excitatory Neuron",
      grepl("GABAergic neuron", sctype_label, ignore.case = TRUE) ~ "GABAergic Neuron",
      grepl("astrocyte", sctype_label, ignore.case = TRUE) ~ "Astrocyte",
      grepl("oligodendro", sctype_label, ignore.case = TRUE) ~ "Oligodendrocyte",
      TRUE ~ sctype_label # If it doesn't match keep the original name
    ))
  
  Idents(tunga_frog_merged) = "broad_annotation"
  
  
  saveRDS(tunga_frog_merged, "tunga_frog_annotated.rds")
}


#DimPlot(tunga_frog_merged, reduction = "umap_by_pca", label = TRUE, repel = TRUE, label.size = 4) +  ggtitle("Tungara frog brain")




# CHECKPOINT 6: Differential expression by sound
# ==============================================================================


#DefaultAssay(tunga_frog_merged) = "RNA"


file_wc_w  = "tunga_frog_WhineChuck_vs_Whine_DE.csv"
file_wc_gn = "tunga_frog_WhineChuck_vs_GreenNoise_DE.csv"
file_w_gn  = "tunga_frog_Whine_vs_GreenNoise_DE.csv"

if (file.exists(file_wc_w) & file.exists(file_wc_gn) & file.exists(file_w_gn)) {
  
  all_markers_WC_vs_W  = read.csv(file_wc_w)
  all_markers_WC_vs_GN = read.csv(file_wc_gn)
  all_markers_W_vs_GN  = read.csv(file_w_gn)
  
} else {
  
  options(future.globals.maxSize = 20000 * 1024^2)
  plan(multisession, workers = 8)
  
  DefaultAssay(tunga_frog_merged) = "RNA"
  tunga_frog_merged$celltype_stimulus = paste(tunga_frog_merged$broad_annotation, tunga_frog_merged$orig.ident, sep = "_")
  Idents(tunga_frog_merged) = "celltype_stimulus"
  
  cell_types = unique(tunga_frog_merged$broad_annotation)
  
  list_wc_w  = list()
  list_wc_gn = list()
  list_w_gn  = list()
  
  for (cell in cell_types){
    

    ident_wc = paste0(cell, "_Whine-chuck")
    ident_w  = paste0(cell, "_Whine")
    ident_gn = paste0(cell, "_Green-noise")
    
    n_wc = sum(Idents(tunga_frog_merged) == ident_wc)
    n_w  = sum(Idents(tunga_frog_merged) == ident_w)
    n_gn = sum(Idents(tunga_frog_merged) == ident_gn)
    
    if (n_wc >= 3 && n_w >= 3) {
      res = FindMarkers(tunga_frog_merged, ident.1 = ident_wc, ident.2 = ident_w, verbose = FALSE)
      if (nrow(res) > 0) { 
        res$cell_type = cell
        res$gene = rownames(res)
        list_wc_w[[cell]] = res 
        }
    }
    
    if (n_wc >= 3 && n_gn >= 3) {
      res = FindMarkers(tunga_frog_merged, ident.1 = ident_wc, ident.2 = ident_gn, verbose = FALSE)
      if (nrow(res) > 0) { 
        res$cell_type = cell
        res$gene = rownames(res)
        list_wc_gn[[cell]] = res 
        }
    }
    
    if (n_w >= 3 && n_gn >= 3) {
      res = FindMarkers(tunga_frog_merged, ident.1 = ident_w, ident.2 = ident_gn, verbose = FALSE)
      if (nrow(res) > 0) { 
        res$cell_type = cell
        res$gene = rownames(res)
        list_w_gn[[cell]] = res 
        }
    }
  }
  
  plan(sequential)
  
  if (length(list_wc_w) > 0) {
    all_markers_WC_vs_W = do.call(rbind, list_wc_w)
    rownames(all_markers_WC_vs_W) = NULL 
    write.csv(all_markers_WC_vs_W, file_wc_w, row.names = FALSE)
  }
  
  if (length(list_wc_gn) > 0) {
    all_markers_WC_vs_GN = do.call(rbind, list_wc_gn)
    rownames(all_markers_WC_vs_GN) = NULL 
    write.csv(all_markers_WC_vs_GN, file_wc_gn, row.names = FALSE)
  }
  
  if (length(list_w_gn) > 0) {
    all_markers_W_vs_GN = do.call(rbind, list_w_gn)
    rownames(all_markers_W_vs_GN) = NULL 
    write.csv(all_markers_W_vs_GN, file_w_gn, row.names = FALSE)
  }
}




top_genes_w_wc = all_markers_WC_vs_W %>%
  filter(p_val_adj < 0.05) %>%
  filter(abs(avg_log2FC) > 1) %>%
  arrange(desc(avg_log2FC))     

top_genes_w_gn = all_markers_W_vs_GN %>%
  filter(p_val_adj < 0.05) %>%
  filter(abs(avg_log2FC) > 1) %>%
  arrange(desc(avg_log2FC))     

top_genes_wc_gn = all_markers_WC_vs_GN %>%
  filter(p_val_adj < 0.05) %>%
  filter(abs(avg_log2FC) > 1) %>%
  arrange(desc(avg_log2FC))     




# head(top_genes_w_wc, 10) 
# p_val avg_log2FC pct.1 pct.2     p_val_adj                  cell_type         gene
# 1  6.110882e-133   9.594846 0.051 0.000 4.368547e-128           GABAergic Neuron LOC140117055
# 2  1.434765e-218   8.054115 0.345 0.001 1.025685e-213          Excitatory Neuron      CAPRIN2
# 3   0.000000e+00   8.035667 0.313 0.001  0.000000e+00           GABAergic Neuron      CAPRIN2
# 4   2.045961e-20   8.027099 0.068 0.000  1.462616e-15                  Astrocyte LOC140077248
# 5   4.954876e-47   8.008103 0.160 0.001  3.542142e-42                  Astrocyte      CAPRIN2
# 6   4.488368e-08   7.421710 0.051 0.000  3.208644e-03 Capillary endothelial cell      CAPRIN2
# 7   1.030027e-60   7.289675 0.024 0.000  7.363454e-56           GABAergic Neuron        DSCC1
# 8   1.108674e-24   6.665361 0.274 0.003  7.925691e-20                Erythrocyte      CAPRIN2
# 9   2.022219e-83   6.569829 0.033 0.000  1.445644e-78           GABAergic Neuron LOC140119394
# 10  5.635192e-58   6.511788 0.107 0.002  4.028486e-53          Excitatory Neuron        SMYD1
# 
# 

# 
# head(top_genes_w_gn, 10)
# p_val avg_log2FC pct.1 pct.2    p_val_adj                  cell_type         gene
# 1  7.238641e-26  11.842621   0.2 0.000 5.174760e-21 Proliferating myeloid cell       GUCY2F
# 2  1.475848e-13  11.376342   0.2 0.002 1.055054e-08 Proliferating myeloid cell         TPH1
# 3  7.238641e-26  10.792389   0.2 0.000 5.174760e-21 Proliferating myeloid cell       TSPAN1
# 4  2.036967e-09  10.355180   0.2 0.004 1.456187e-04 Proliferating myeloid cell         ASMT
# 5  7.238641e-26   9.685340   0.2 0.000 5.174760e-21 Proliferating myeloid cell        PDE6G
# 6  7.238641e-26   9.569368   0.2 0.000 5.174760e-21 Proliferating myeloid cell LOC140120519
# 7  7.238641e-26   9.569368   0.2 0.000 5.174760e-21 Proliferating myeloid cell      SLCO2A1
# 8  7.238641e-26   9.569368   0.2 0.000 5.174760e-21 Proliferating myeloid cell LOC140128939
# 9  7.238641e-26   9.569368   0.2 0.000 5.174760e-21 Proliferating myeloid cell        ITIH2
# 10 7.238641e-26   9.569368   0.2 0.000 5.174760e-21 Proliferating myeloid cell LOC140104537

# head(top_genes_wc_gn, 10)
# p_val avg_log2FC pct.1 pct.2     p_val_adj                  cell_type         gene
# 1   1.474840e-29   9.536801 0.051 0.000  1.054333e-24 Capillary endothelial cell      CAPRIN2
# 2   4.868127e-41   9.422197 0.161 0.000  3.480127e-36                  Astrocyte LOC140126440
# 3   2.080354e-68   8.542374 0.226 0.002  1.487203e-63            Oligodendrocyte LOC140116741
# 4  1.475701e-194   7.734227 0.345 0.002 1.054949e-189          Excitatory Neuron      CAPRIN2
# 5   8.565863e-08   7.664836 0.011 0.000  6.123564e-03 Capillary endothelial cell LOC140121587
# 6   2.758523e-07   7.646857 0.078 0.000  1.972013e-02    Antigen-presenting cell LOC140069429
# 7   8.565863e-08   7.518484 0.011 0.000  6.123564e-03 Capillary endothelial cell LOC140069072
# 8   0.000000e+00   7.456371 0.313 0.002  0.000000e+00           GABAergic Neuron      CAPRIN2
# 9   7.726548e-28   7.345055 0.095 0.001  5.523554e-23            Oligodendrocyte      CAPRIN2
# 10  6.293120e-24   6.949246 0.010 0.000  4.498825e-19           GABAergic Neuron LOC140125685


 VlnPlot(tunga_frog_merged, features = "CAPRIN2", group.by = "orig.ident")
 # VlnPlot(tunga_frog_merged, features = "LOC140117055", group.by = "orig.ident")
 # VlnPlot(tunga_frog_merged, features = "LOC140077248", group.by = "orig.ident")
 # VlnPlot(tunga_frog_merged, features = "DSCC1", group.by = "orig.ident")
 # VlnPlot(tunga_frog_merged, features = "LOC140119394", group.by = "orig.ident")
 # VlnPlot(tunga_frog_merged, features = "SMYD1", group.by = "orig.ident")
 
 VlnPlot(tunga_frog_merged, features = "GUCY2F", group.by = "orig.ident")
 
 VlnPlot(tunga_frog_merged, features = "LOC140126440", group.by = "orig.ident")
 
 
 # CHECKPOINT 7: Map LOC genes to Xenopus orthologs
 # ==============================================================================

  
 orthos = read.csv("xenopus_gene_to_tungara_gene.csv")
 map_loc_genes = function(de_results, orthos_df) {

   unique_orthos = orthos_df %>% distinct(tungara_gene, xenopus_gene, .keep_all = TRUE)
   
   de_results %>% left_join(unique_orthos, by = c("gene" = "tungara_gene")) %>%
     mutate(translated_gene = coalesce(xenopus_gene, gene)) %>%
     relocate(translated_gene, .after = gene) 
 }
 
 all_markers_WC_vs_W  = map_loc_genes(all_markers_WC_vs_W, orthos)
 all_markers_WC_vs_GN = map_loc_genes(all_markers_WC_vs_GN, orthos)
 all_markers_W_vs_GN  = map_loc_genes(all_markers_W_vs_GN, orthos)
 
 top_genes_w_wc = all_markers_WC_vs_W %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 1) %>% arrange(desc(avg_log2FC))
 top_genes_w_gn = all_markers_W_vs_GN %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 1) %>% arrange(desc(avg_log2FC))
 top_genes_wc_gn = all_markers_WC_vs_GN %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 1) %>% arrange(desc(avg_log2FC))
 
 
 head(top_genes_w_wc, 10) 

 
 
 
 # CHECKPOINT 8: hdWGCNA
 # ==============================================================================
 
#install.packages("BiocManager")
#BiocManager::install(c("harmony", "UCell", "GeneOverlap"))
#devtools::install_github("immunogenomics/harmony")
#devtools::install_github('smorabit/hdWGCNA', ref='dev')

 
 if (file.exists("tunga_frog_annotated_wgcna.rds")) {
   
   tunga_frog_merged = readRDS("tunga_frog_annotated_wgcna.rds")
   
 } else { 
  
  tunga_frog_merged = SetupForWGCNA(
     tunga_frog_merged,
     gene_select = "fraction", # Select genes expressed in a certain fraction of cells
     fraction = 0.05,          # In at least 5% of cells
     wgcna_name = "Tungara_WGCNA"
   )
   
  tunga_frog_merged = MetacellsByGroups(
     seurat_obj = tunga_frog_merged,
     group.by = c("broad_annotation", "orig.ident"), 
     k = 25, # k-nearest neighbors 
     max_shared = 10, 
     ident.group = "broad_annotation"
   )
   
  tunga_frog_merged = NormalizeMetacells(tunga_frog_merged)
  
  metacell_obj = GetMetacellObject(tunga_frog_merged)
  
  valid_cell_types = as.character(unique(metacell_obj$broad_annotation))
  
  # Expression matrix for WGCNA 
  tunga_frog_merged = SetDatExpr(
     tunga_frog_merged,
     group_name = valid_cell_types, 
     group.by = "broad_annotation",
     assay = "RNA",
     slot = "data"
   )
   
  # Soft-thresholding powers to have scale-free topology
  tunga_frog_merged = TestSoftPowers(
     tunga_frog_merged,
     networkType = 'signed' # Signed networks
   )
  
  PlotSoftPowers(tunga_frog_merged) 
  
  # Construct the co-expression network 
  tunga_frog_merged = ConstructNetwork(
     tunga_frog_merged, 
     soft_power = 10, # Power from TestSoftPowers
     setDatExpr = FALSE,
     tom_name = 'Tungara_TOM' 
   )
   
   # Compute module eigengenes and module connectivity
  tunga_frog_merged = ModuleEigengenes(tunga_frog_merged)
  tunga_frog_merged = ModuleConnectivity(tunga_frog_merged)
  
  eigengenes = GetMEs(tunga_frog_merged)
  tunga_frog_merged = AddMetaData(tunga_frog_merged, metadata = eigengenes)
  
  
  saveRDS(tunga_frog_merged, "tunga_frog_annotated_wgcna.rds")
 }
 
 
PlotDendrogram(tunga_frog_merged, main = 'Tungara Frog Brain Gene Network')

ModuleFeaturePlot(
  tunga_frog_merged,
  features = 'hMEs', # Module Eigengenes
  order = TRUE,       # cells with the highest expression to the front
  reduction = 'umap_by_pca'
)

head(GetModules(tunga_frog_merged))

options(future.globals.maxSize = 20000 * 1024^2)

#graphics.off()

HubGeneNetworkPlot(
  tunga_frog_merged,
  n_hubs = 25,       # 25 most connected genes
  n_other = 5,       # 5 loosely connected genes for context
  edge_prop = 0.75,  # thickness of the lines
  mods = "brown" # Specific module color
)

tunga_neurons = subset(tunga_frog_merged, idents = c("Excitatory Neuron", "GABAergic Neuron"))

VlnPlot(
  tunga_neurons,
  features = "brown", 
  group.by = "orig.ident", 
  pt.size = 0
) + ggtitle("Brown: Neurons Only")

