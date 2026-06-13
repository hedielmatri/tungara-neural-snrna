setwd("/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/Hedi_R")

library(Seurat)
library(dplyr)
library(future)
set.seed(123)

# CHECKPOINT 1: Tungara Data Processing
# ==============================================================================
if (file.exists("/tunga_frog_clustered.rds")) {
  tunga_frog_merged = readRDS("/tunga_frog_clustered.rds")
} else {
  filtered_counts_gn = Read10X_h5("stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/TF_GN_outs/outs/filtered_feature_bc_matrix.h5")
  filtered_counts_wc = Read10X_h5("stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/TF_GN_outs/outs/TF_WC_outs/outs/filtered_feature_bc_matrix.h5")
  filtered_counts_w  = Read10X_h5("stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/TF_GN_outs/outs/TF_W_outs/outs/filtered_feature_bc_matrix.h5")
  
  seurat_gn = CreateSeuratObject(counts = filtered_counts_gn, project = "Green-noise")
  seurat_wc = CreateSeuratObject(counts = filtered_counts_wc, project = "Whine-chuck")
  seurat_w  = CreateSeuratObject(counts = filtered_counts_w, project = "Whine")
  
  tunga_frog_merged = merge(x = seurat_gn, y = c(seurat_wc, seurat_w), add.cell.ids = c("GN", "WC", "W"), project = "Tunga_Frog")
  tunga_frog_merged = subset(tunga_frog_merged, subset = nFeature_RNA > 500 & nFeature_RNA < 6000)
  
  tunga_frog_merged = NormalizeData(tunga_frog_merged) 
  tunga_frog_merged = FindVariableFeatures(tunga_frog_merged) 
  tunga_frog_merged = ScaleData(tunga_frog_merged) 
  tunga_frog_merged = RunPCA(tunga_frog_merged)
  
  
  tunga_frog_merged = IntegrateLayers(
    object = tunga_frog_merged, method = CCAIntegration, 
    orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE
  )
  
  tunga_frog_merged = FindNeighbors(tunga_frog_merged, dims = 1:20, reduction = "pca") 
  tunga_frog_merged = FindClusters(tunga_frog_merged, resolution = 0.2) 
  tunga_frog_merged = RunUMAP(tunga_frog_merged, dims = 1:20, reduction = "pca", reduction.name = "umap_by_pca")
  tunga_frog_merged = JoinLayers(tunga_frog_merged)
  
  saveRDS(tunga_frog_merged, "/tunga_frog_clustered.rds") 
}


# CHECKPOINT 2: FindAllMarkers
# ==============================================================================
if (file.exists("/tunga_frog_markers.rds")) {
  markers = readRDS("/tunga_frog_markers.rds")
} else {
  plan(multisession, workers = 22)
  markers = FindAllMarkers(
    tunga_frog_merged, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25,
    max.cells.per.ident = 500
  )
  plan(sequential)
  saveRDS(markers, "/tunga_frog_markers.rds") 
}

# CHECKPOINT 3: Xenopus reference
# ==============================================================================
if (file.exists("/xeno_brain_clean.rds")) {
  xeno_brain_clean = readRDS("/xeno_brain_clean.rds")
  xeno_meta = read.csv("/Brain_cell_info.csv", row.names = 1) 
} else {
  brain1 = read.table("GSM6214268_Xenopus_brain_COL65_dge.txt", header = TRUE, row.names = 1)
  brain2 = read.table("GSM6214269_Xenopus_brain_COL66_dge.txt", header = TRUE, row.names = 1)
  brain3 = read.table("GSM6214270_Xenopus_brain_COL124_dge.txt", header = TRUE, row.names = 1)
  brain4 = read.table("GSM6214271_Xenopus_brain_COL127_dge.txt", header = TRUE, row.names = 1)
  
  seu_b1 = CreateSeuratObject(counts = brain1, project = "XenoBrain")
  seu_b2 = CreateSeuratObject(counts = brain2, project = "XenoBrain")
  seu_b3 = CreateSeuratObject(counts = brain3, project = "XenoBrain")
  seu_b4 = CreateSeuratObject(counts = brain4, project = "XenoBrain")
  xeno_brain = merge(seu_b1, y = c(seu_b2, seu_b3, seu_b4), add.cell.ids = c("b1", "b2", "b3", "b4"))
  
  xeno_meta = read.csv("/Brain_cell_info.csv", row.names = 1)
  
  # --- THE ROBUST METADATA FIX ---
  seurat_batch_raw = sub("_.*", "", colnames(xeno_brain))  # Extracts "b1", "b2"
  seurat_dna = sub(".*_", "", colnames(xeno_brain))        # Extracts DNA barcode
  seurat_batch = sub("b", "Brain", seurat_batch_raw)       # Converts "b1" to "Brain1"
  seurat_keys = paste0(seurat_batch, "_", seurat_dna)      # Creates "Brain1_AACCT..."
  
  meta_batch = sub("_.*", "", rownames(xeno_meta))         # Extracts "Brain1"
  meta_dna = sub(".*\\.", "", rownames(xeno_meta))         # Extracts DNA barcode
  meta_keys = paste0(meta_batch, "_", meta_dna)            # Creates "Brain1_AACCT..."
  
  # Match and assign
  match_indices = match(seurat_keys, meta_keys)
  xeno_brain$cluster = xeno_meta$cluster[match_indices]
  xeno_brain = subset(xeno_brain, subset = !is.na(cluster))
  # -------------------------------
  
  # Clean the Reference (Pituitary removal)
  cluster_to_celltype = unique(xeno_meta[, c("cluster", "celltype")])
  name_map = setNames(cluster_to_celltype$celltype, as.character(cluster_to_celltype$cluster))
  xeno_brain$celltype = unname(name_map[as.character(xeno_brain$cluster)])
  
  pituitary_keywords = "Gonadotroph|Thyrotroph|Endocrine|Growth hormone|Prolactin|Melanotrope"
  xeno_brain_clean = subset(xeno_brain, subset = !grepl(pituitary_keywords, celltype))
  
  saveRDS(xeno_brain_clean, "/xeno_brain_clean.rds") # 
}

# 4. TRANSLATE THE GENOME
# ==============================================================================
tung_counts = LayerData(tunga_frog_merged, assay = "RNA", layer = "counts")
rownames(tung_counts) = tolower(rownames(tung_counts))

xeno_counts = LayerData(xeno_brain_clean, assay = "RNA", layer = "counts")
current_xeno_genes = tolower(sub("\\..*", "", rownames(xeno_counts)))

orthos = read.csv("xenopus_gene_to_tungara_gene.csv")
orthos$xenopus_gene = tolower(orthos$xenopus_gene)
orthos$tungara_gene = tolower(orthos$tungara_gene)

patched_genes = current_xeno_genes
match_idx = match(current_xeno_genes, orthos$xenopus_gene)
has_translation = !is.na(match_idx)
patched_genes[has_translation] = orthos$tungara_gene[match_idx[has_translation]]

valid_idx = !duplicated(patched_genes) & !is.na(patched_genes)
xeno_counts_translated = xeno_counts[valid_idx, ]
rownames(xeno_counts_translated) = patched_genes[valid_idx]

common_genes = intersect(rownames(tung_counts), rownames(xeno_counts_translated))

tung_counts_common = tung_counts[common_genes, ]
tunga_frog_merged[["XenoMapped"]] = CreateAssayObject(counts = tung_counts_common)
DefaultAssay(tunga_frog_merged) = "XenoMapped"
tunga_frog_merged = NormalizeData(tunga_frog_merged)

xeno_counts_common = xeno_counts_translated[common_genes, ]
xeno_brain_mapped = CreateSeuratObject(counts = xeno_counts_common, meta.data = xeno_brain_clean@meta.data)
xeno_brain_mapped = NormalizeData(xeno_brain_mapped)

# CHECKPOINT 4: Label transfer
# ==============================================================================
if (file.exists("/tunga_frog_integrated.rds")) {
  tunga_frog_merged = readRDS("/tunga_frog_integrated.rds")
} else {
  anchors = FindTransferAnchors(
    reference = xeno_brain_mapped,
    query = tunga_frog_merged,
    features = common_genes,
    reference.assay = "RNA",
    query.assay = "XenoMapped",
    reduction = "cca"
  )
  
  predictions = TransferData(
    anchorset = anchors,
    refdata = as.character(xeno_brain_mapped$cluster),
    weight.reduction = tunga_frog_merged[["pca"]], 
    dims = 1:20
  )
  
  tunga_frog_merged = AddMetaData(tunga_frog_merged, metadata = predictions$predicted.id, col.name = "Xeno_Brain_Clusters")
  
  cluster_to_celltype = unique(xeno_meta[, c("cluster", "celltype")])
  name_map = setNames(cluster_to_celltype$celltype, as.character(cluster_to_celltype$cluster))
  tunga_frog_merged$Xeno_CellTypes = unname(name_map[as.character(tunga_frog_merged$Xeno_Brain_Clusters)])
  
  majority_vote = tunga_frog_merged@meta.data %>% 
    count(seurat_clusters, Xeno_CellTypes) %>% 
    group_by(seurat_clusters) %>% 
    slice_max(n, n = 1, with_ties = FALSE)
  
  final_mapping = setNames(majority_vote$Xeno_CellTypes, majority_vote$seurat_clusters)
  Idents(tunga_frog_merged) = "seurat_clusters"
  tunga_frog_merged = RenameIdents(tunga_frog_merged, final_mapping)
  tunga_frog_merged$Final_Annotations = Idents(tunga_frog_merged)
  
  saveRDS(tunga_frog_merged, "/tunga_frog_integrated.rds") 
}


# 5. FINAL PLOT
# ==============================================================================
DimPlot(tunga_frog_merged, reduction = "umap_by_pca", label = TRUE, label.size = 3) + ggtitle("Tungara frog brain (Final Clean Annotation)")
