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

# CHECKPOINT 9 : Annotation by Blast
-
if (file.exists("tunga_frog_annotated_blast.rds")) {
  tunga_frog_merged_blast = readRDS("tunga_frog_annotated_blast.rds")
} else {
  tunga_frog_merged_blast = readRDS("tunga_frog_clustered.rds")
  
  brain = read_excel("xenopus_brain_markers.xlsx", sheet = "Brain")
  brain = brain %>% tidyr::fill(`Cell type`) 
  

  # LOAD BLAST MATCHES
  orthos = read.csv("xenopus_to_tungara_RBH.csv")
  orthos$xenopus_protein = trimws(gsub("\\..*", "", orthos$xenopus_protein))
  orthos$tungara_protein = trimws(gsub("\\..*", "", orthos$tungara_protein))
  print(paste(nrow(orthos), "protein matches"))
  
  # LOAD THE GOLDEN TROPICALIS DICT FOUND
  xeno_id_to_symbol = read.csv("xenopus_id_to_symbol_tropicalis.csv") %>% 
    dplyr::rename(xenopus_protein = protein_id, xenopus_symbol = gene_symbol) %>%
    dplyr::filter(grepl("^[XN]P_", xenopus_protein)) 
  
  xeno_id_to_symbol$xenopus_protein = trimws(gsub("\\..*", "", xeno_id_to_symbol$xenopus_protein))
  
  # LOAD THE REFSEQ DICTIONARY
  tunga_id_to_symbol = read.csv("tungara_id_to_symbol_FIXED.csv") %>% 
    dplyr::rename(tungara_protein = protein_id, tungara_symbol = gene_symbol)  
  
  tunga_id_to_symbol$tungara_protein = trimws(gsub("\\..*", "", tunga_id_to_symbol$tungara_protein))
  
  # Build the Master Translation Table
  translation_table = orthos %>%
    dplyr::inner_join(xeno_id_to_symbol, by = "xenopus_protein") %>%
    dplyr::inner_join(tunga_id_to_symbol, by = "tungara_protein") %>%
    dplyr::mutate(
      xenopus_clean = toupper(xenopus_symbol),
      tungara_symbol = toupper(tungara_symbol)
    ) %>%
    dplyr::distinct(xenopus_clean, tungara_symbol, .keep_all = TRUE)
  
  # write_xlsx(translation_table, "translation_table_xenopus_to_tungara.xlsx")

  print(paste("TRANSLATION TABLE:", nrow(translation_table), "matchs"))
  seurat_genes = toupper(rownames(tunga_frog_merged_blast))
  
  # Apply the translation to the Excel markers
  brain_mapped = brain %>%
    dplyr::mutate(xenopus_clean = toupper(gsub("\\..*", "", gene))) %>%
    dplyr::inner_join(translation_table, by = "xenopus_clean") %>%
    dplyr::filter(tungara_symbol %in% seurat_genes) %>%
    dplyr::distinct(`Cell type`, gene, tungara_symbol, .keep_all = TRUE)
  
  print(paste("Number of valid brain markers mapped to Seurat:", nrow(brain_mapped)))
  
  # Filter strong markers
  brain_mapped_strong = brain_mapped %>%
    dplyr::group_by(`Cell type`) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) %>%
    dplyr::slice_head(n = 50)
  
  # Build ScType dictionary
  sctype_db = brain_mapped_strong %>%
    dplyr::group_by(`Cell type`) %>%
    dplyr::summarize(
      geneSymbolmore1 = paste(unique(tungara_symbol), collapse = ","),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      tissueType = "Brain",
      geneSymbolmore2 = ""
    ) %>%
    dplyr::rename(cellName = `Cell type`) %>%
    dplyr::select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2)
  
  openxlsx::write.xlsx(sctype_db, "tungara_brain_sctype_markers_clean.xlsx", rowNames = FALSE)
  
  print("Running ScType annotations")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  gs_list = gene_sets_prepare("tungara_brain_sctype_markers_clean.xlsx", "Brain")
  
  # Raw compressed data
  expr_matrix = GetAssayData(tunga_frog_merged_blast, assay = "RNA", slot = "data")
  
  # ONLY marker genes ScType need
  marker_genes = unique(c(unlist(gs_list$gs_positive), unlist(gs_list$gs_negative)))
  valid_markers = marker_genes[marker_genes %in% rownames(expr_matrix)]
  
  # Matrix from 30 000 genes down to just the valid markers
  expr_matrix_small = expr_matrix[valid_markers, ]
  
  print(paste("Matrix of", nrow(expr_matrix_small), "genes"))
  
  # ScType
  es.max = sctype_score(
    scRNAseqData = expr_matrix_small,
    scaled = FALSE, 
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  
  # Apply labels
  sctype_top = do.call("rbind", lapply(unique(tunga_frog_merged_blast$seurat_clusters), function(cl) {
    cells = rownames(tunga_frog_merged_blast@meta.data)[tunga_frog_merged_blast$seurat_clusters == cl]
    scores = sort(rowSums(es.max[, cells, drop = FALSE]), decreasing = TRUE)
    
    data.frame(
      cluster = cl,
      cell_type = names(scores)[1],
      score = scores[1]
    )
  }))
  
  tunga_frog_merged_blast$sctype_label = sctype_top$cell_type[
    match(tunga_frog_merged_blast$seurat_clusters, sctype_top$cluster)
  ]
  
  tunga_frog_merged_blast@meta.data = tunga_frog_merged_blast@meta.data %>%
    dplyr::mutate(broad_annotation = dplyr::case_when(
      grepl("Excitatory neuron", sctype_label, ignore.case = TRUE) ~ "Excitatory Neuron",
      grepl("GABAergic neuron", sctype_label, ignore.case = TRUE) ~ "GABAergic Neuron",
      grepl("astrocyte", sctype_label, ignore.case = TRUE) ~ "Astrocyte",
      grepl("oligodendro", sctype_label, ignore.case = TRUE) ~ "Oligodendrocyte",
      TRUE ~ sctype_label 
    ))
  
  Idents(tunga_frog_merged_blast) = "broad_annotation"
  
  saveRDS(tunga_frog_merged_blast, "tunga_frog_annotated_blast.rds")
}



DimPlot(tunga_frog_merged_blast, reduction = "umap_by_pca", label = TRUE, repel = TRUE, label.size = 4) +  ggtitle("Tungara frog brain blast")

table(tunga_frog_merged_blast$broad_annotation)

print(sum(table(tunga_frog_merged_blast$broad_annotation)))




#Idents(tunga_frog_merged_blast) = "broad_annotation"
#markers = FindAllMarkers(tunga_frog_merged_blast, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1.0)
#top1 = markers %>% group_by(cluster) %>% slice_head(n = 5)

#DotPlot(tunga_frog_merged_blast, features = unique(top1$gene)) + RotatedAxis() + ggtitle("Cell type markers")







# CHECKPOINT 10: GABA SUBCLUSTER ANNOTATION

if (file.exists("tunga_gaba_subclusters_annotated.rds")) {
  
  gaba_subset = readRDS("tunga_gaba_subclusters_annotated.rds")
  
  
} else {
  
  
  
  # ONLY the GABAergic cluster
  gaba_subset = subset(tunga_frog_merged_blast, idents = "GABAergic Neuron")
  
  gaba_subset = FindVariableFeatures(gaba_subset, selection.method = "vst", nfeatures = 2000)
  
  # only the top 2000 genes
  gaba_subset = ScaleData(gaba_subset)
  
  
  gaba_subset = RunPCA(gaba_subset, features = VariableFeatures(object = gaba_subset), verbose = FALSE)
  gaba_subset = FindNeighbors(gaba_subset, dims = 1:20, verbose = FALSE)
  
  # resolution of 0.4 to force it to break down
  gaba_subset = FindClusters(gaba_subset, resolution = 0.4, verbose = FALSE)
  gaba_subset = RunUMAP(gaba_subset, dims = 1:20, verbose = FALSE)
  
  #saveRDS(gaba_subset, "tunga_gaba_subclusters.rds")
  #gaba_subset = readRDS("tunga_gaba_subclusters.rds")
  
  #DimPlot(gaba_subset, reduction = "umap", label = TRUE, label.size = 5) + ggtitle("GABAergic broken down")
  
  
  gaba_markers = FindAllMarkers(gaba_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2)
  
  # 2. Build a reference dictionary from the BLAST-mapped Xenopus data
  blast_sctype_db = brain_mapped %>%
    dplyr::group_by(`Cell type`) %>%
    dplyr::slice_max(n = 50, order_by = avg_log2FC) %>%
    dplyr::summarize(
      geneSymbolmore1 = paste(unique(tungara_symbol), collapse = ","),
      .groups = "drop"
    ) %>%
    dplyr::mutate(tissueType = "Brain", geneSymbolmore2 = "") %>%
    dplyr::rename(cellName = `Cell type`) %>%
    dplyr::select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2)
  
  openxlsx::write.xlsx(blast_sctype_db, "gaba_subclusters_ref_markers.xlsx", rowNames = FALSE)
  
  # Score the sub-clusters against reference
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  gs_list_ref = gene_sets_prepare("gaba_subclusters_ref_markers.xlsx", "Brain")
  
  expr_matrix_gaba = GetAssayData(gaba_subset, assay = "RNA", slot = "data")
  
  es.max.ref = sctype_score(
    scRNAseqData = expr_matrix_gaba[rownames(expr_matrix_gaba) %in% unlist(gs_list_ref), ],
    scaled = FALSE, 
    gs = gs_list_ref$gs_positive,
    gs2 = gs_list_ref$gs_negative
  )
  
  # Assign labels based 
  gaba_labels = do.call("rbind", lapply(unique(gaba_subset$seurat_clusters), function(cl) {
    cells = rownames(gaba_subset@meta.data)[gaba_subset$seurat_clusters == cl]
    scores = sort(rowSums(es.max.ref[, cells, drop = FALSE]), decreasing = TRUE)
    data.frame(cluster = cl, cell_type = names(scores)[1], score = scores[1])
  }))
  
  # Apply labels
  gaba_subset$gaba_ref_annotation = gaba_labels$cell_type[match(gaba_subset$seurat_clusters, gaba_labels$cluster)]
  Idents(gaba_subset) = "gaba_ref_annotation"
  
  gaba_subset@meta.data = gaba_subset@meta.data %>%
    dplyr::mutate(gaba_broad_annotation = dplyr::case_when(
      grepl("Excitatory neuron", gaba_ref_annotation, ignore.case = TRUE) ~ "Excitatory Neuron",
      grepl("GABAergic neuron", gaba_ref_annotation, ignore.case = TRUE) ~ "GABAergic Neuron",
      grepl("Growth hormone", gaba_ref_annotation, ignore.case = TRUE) ~ "Growth hormone",
      TRUE ~ gaba_ref_annotation 
    ))
  
  Idents(gaba_subset) = "gaba_ref_annotation"
  
  saveRDS(gaba_subset, "tunga_gaba_subclusters_annotated.rds")
  
}


DimPlot(gaba_subset, reduction = "umap", label = TRUE, repel = TRUE) +  ggtitle("GABAergic subsets annotated")

# Check the table
print(table(gaba_subset$gaba_ref_annotation))



# CHECKPOINT 11 : FINAL OBJECT
# ==============================================================================

if (file.exists("tunga_frog_blast_with_gaba_subtypes.rds")) {
  
  tunga_frog_merged_blast = readRDS( "tunga_frog_blast_with_gaba_subtypes.rds")
  
  
} else {
  
  
  
  tunga_frog_merged_blast$refined_annotation = as.character(tunga_frog_merged_blast$sctype_label)
  
  Idents(gaba_subset) = "gaba_ref_annotation"
  #Idents(gaba_subset) = "gaba_broad_annotation"
  
  refined_labels = Idents(gaba_subset)
  matching_cells = names(refined_labels)
  
  tunga_frog_merged_blast$refined_annotation[match(matching_cells, rownames(tunga_frog_merged_blast@meta.data))] = as.character(refined_labels)
  
  Idents(tunga_frog_merged_blast) = "refined_annotation"
  
  
  saveRDS(tunga_frog_merged_blast, "tunga_frog_blast_with_gaba_subtypes.rds")
}



DimPlot(tunga_frog_merged_blast, reduction = "umap_by_pca", label = TRUE, repel = TRUE, label.size = 3) + 
  ggtitle("Brain by Blast with GABAergic subtypes")


if (file.exists("tunga_frog_annotated_blast_wgcna.rds")) {
  
  tunga_frog_merged_blast = readRDS("tunga_frog_annotated_blast_wgcna.rds")
  
} else { 
  
  
  tunga_frog_merged_blast@meta.data = tunga_frog_merged_blast@meta.data %>%
    dplyr::mutate(blast_gabbasub_broad_annotation = dplyr::case_when(
      grepl("Excitatory neuron", refined_annotation, ignore.case = TRUE) ~ "Excitatory Neuron",
      grepl("GABAergic neuron", refined_annotation, ignore.case = TRUE) ~ "GABAergic Neuron",
      
      grepl("astrocyte", refined_annotation, ignore.case = TRUE) ~ "Astrocyte",
      grepl("oligodendro", refined_annotation, ignore.case = TRUE) ~ "Oligodendrocyte",
      
      grepl("endocrine|gonadotroph|thyrotroph|melanotrope|growth hormone|prolactin", refined_annotation, ignore.case = TRUE) ~ "Endocrine/Pituitary Cell",
      
      grepl("endothelial", refined_annotation, ignore.case = TRUE) ~ "Endothelial Cell",
      grepl("antigen-presenting", refined_annotation, ignore.case = TRUE) ~ "Immune Cell",
      
      TRUE ~ refined_annotation 
    ))  
  
  Idents(tunga_frog_merged_blast) = "blast_gabbasub_broad_annotation"
  
  #DimPlot(tunga_frog_merged_blast, reduction = "umap_by_pca", label = TRUE, repel = TRUE, label.size = 3) + ggtitle("Brain by Blast with GABAergic subtypes")
  
  
  tunga_frog_merged_blast = SetupForWGCNA(
    tunga_frog_merged_blast,
    gene_select = "fraction", # Select genes expressed in a certain fraction of cells
    fraction = 0.05,          # In at least 5% of cells
    wgcna_name = "Tungara_WGCNA"
  )
  
  tunga_frog_merged_blast = MetacellsByGroups(
    seurat_obj = tunga_frog_merged_blast,
    group.by = c("blast_gabbasub_broad_annotation", "orig.ident"), 
    k = 25, 
    max_shared = 10, 
    ident.group = "blast_gabbasub_broad_annotation" 
  )
  
  tunga_frog_merged_blast = NormalizeMetacells(tunga_frog_merged_blast)
  
  # Expression for only GABAergic
  tunga_frog_merged_blast = SetDatExpr(
    tunga_frog_merged_blast,
    group_name = "GABAergic Neuron", 
    group.by = "blast_gabbasub_broad_annotation", 
    assay = "RNA",
    slot = "data"
  )  
  
  
  # metacell_obj = GetMetacellObject(tunga_frog_merged_blast)
  # valid_cell_types = as.character(unique(metacell_obj$refined_annotation))
  # 
  # # Expression matrix for WGCNA whole brain
  # tunga_frog_merged_blast = SetDatExpr(
  #   tunga_frog_merged_blast,
  #   group_name = valid_cell_types, 
  #   group.by = "refined_annotation",
  #   assay = "RNA",
  #   slot = "data"
  # )

  tunga_frog_merged_blast = TestSoftPowers(
    tunga_frog_merged_blast,
    networkType = 'signed' # Signed networks
  )
  
  PlotSoftPowers(tunga_frog_merged_blast) 
  
  # Construct co-expression network 
  tunga_frog_merged_blast = ConstructNetwork(
    tunga_frog_merged_blast, 
    soft_power = 7, # Power from TestSoftPowers
    setDatExpr = FALSE,
    tom_name = 'Tungara_TOM',
    overwrite_tom = TRUE, # replacing old
    
    deepSplit = 4,         # split into smaller modules
    minModuleSize = 50,    # Prevents  from being too small
    mergeCutHeight = 0.05   # Merges modules that are >80% correlated
  )
  
  # Compute module eigengenes and module connectivity
  tunga_frog_merged_blast = ModuleEigengenes(tunga_frog_merged_blast)
  tunga_frog_merged_blast = ModuleConnectivity(tunga_frog_merged_blast)
  
  eigengenes = GetMEs(tunga_frog_merged_blast)
  tunga_frog_merged_blast = AddMetaData(tunga_frog_merged_blast, metadata = eigengenes)
  
  
  saveRDS(tunga_frog_merged_blast, "tunga_frog_annotated_blast_wgcna.rds")
}
