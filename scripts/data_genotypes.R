setwd("/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/Hedi_R")

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Matrix)

dir = "/stor/work/Hofmann/All_projects/E_pustulosus_Mate-Choice_Walkowski/JA26095_sequences/souporcell/"

tunga = readRDS("tunga_frog_annotated_blast_wgcna.rds")
tunga = subset(tunga, subset = orig.ident %in% c("Whine", "Whine-chuck"))
tunga$sound = tunga$orig.ident

sc_w = read.table(paste0(dir, "souporcell_W_k2/clusters.tsv"), header = TRUE) %>%
  filter(status == "singlet") %>% mutate(frog = paste0("W", assignment))

sc_wc = read.table(paste0(dir, "souporcell_WC_k2/clusters.tsv"), header = TRUE) %>%
  filter(status == "singlet") %>% mutate(frog = paste0("WC", assignment))

bc = colnames(tunga)
strip_bc = sub("^[A-Za-z]+_", "", bc)
lookup = rep(NA, length(bc))

w_idx = tunga$orig.ident == "Whine"
wc_idx = tunga$orig.ident == "Whine-chuck"

lookup[w_idx] = sc_w$frog[match(strip_bc[w_idx], sc_w$barcode)]
lookup[wc_idx] = sc_wc$frog[match(strip_bc[wc_idx], sc_wc$barcode)]

tunga$frog = lookup
tunga = subset(tunga, !is.na(frog))

counts = GetAssayData(tunga, layer = "counts")
frogs = unique(tunga$frog)

bulk = sapply(frogs, function(f) rowSums(counts[, tunga$frog == f]))
colnames(bulk) = frogs

cpm = sweep(bulk, 2, colSums(bulk), "/") * 1e6
cpm = cpm[rowSums(bulk) >= 20, ]
cpm = cpm[!grepl("^LOC", rownames(cpm)), ]

w_cols = grep("^W", colnames(cpm))
wc_cols = grep("^WC", colnames(cpm))

avg_w = rowMeans(cpm[, w_cols])
avg_wc = rowMeans(cpm[, wc_cols])
fc = log2((avg_wc + 1) / (avg_w + 1))

# SPLIT TOP 40 INTO FIRST 20 AND SECOND 20
sorted_genes = names(sort(abs(fc), decreasing = TRUE))
top20_1 = sorted_genes[1:20]
top20_2 = sorted_genes[21:40]

# FIRST 20 GENES
df1a = data.frame(
  gene = factor(top20_1, levels = top20_1[order(abs(fc[top20_1]))]),
  fc = fc[top20_1],
  direction = ifelse(fc[top20_1] > 0, "Up in WC", "Down in WC")
)
fig1a = ggplot(df1a, aes(x = fc, y = gene, fill = direction)) + geom_col()
print(fig1a)

# SECOND 20 GENES
df1b = data.frame(
  gene = factor(top20_2, levels = top20_2[order(abs(fc[top20_2]))]),
  fc = fc[top20_2],
  direction = ifelse(fc[top20_2] > 0, "Up in WC", "Down in WC")
)
fig1b = ggplot(df1b, aes(x = fc, y = gene, fill = direction)) + geom_col()
print(fig1b)

# FIRST 20 GENES BY FROG
df2a = as.data.frame(cpm[top20_1, ]) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "frog", values_to = "cpm") %>%
  mutate(gene = factor(gene, levels = top20_1))

fig2a = ggplot(df2a, aes(x = frog, y = cpm, group = frog)) +
  geom_point() + geom_line() +
  facet_wrap(~ gene, scales = "free_y")
print(fig2a)

# SECOND 20 GENES BY FROG
df2b = as.data.frame(cpm[top20_2, ]) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "frog", values_to = "cpm") %>%
  mutate(gene = factor(gene, levels = top20_2))

fig2b = ggplot(df2b, aes(x = frog, y = cpm, group = frog)) +
  geom_point() + geom_line() +
  facet_wrap(~ gene, scales = "free_y")
print(fig2b)

# CONSISTENT GENES FROM THE TOP 20 LIST
final_genes = c("CAPRIN2", "CPZ", "ENDOV", "PLG", "SMYD1", "LCA5L", "ABCB10", "SH3YL1")
final_genes = final_genes[final_genes %in% rownames(cpm)]

final_df = cpm[final_genes, ]
print(final_df)

# FIGURE 3
tunga$cell_group = ifelse(tunga$broad_annotation %in% c("Excitatory Neuron", "GABAergic Neuron"), "Neuron", "Other")

f3_list = lapply(final_genes, function(g) {
  d = data.frame(expr = FetchData(tunga, vars = g)[,1] > 0,
                 group = tunga$cell_group, 
                 sound = tunga$sound)
  res = d %>% group_by(group, sound) %>% summarize(pct = mean(expr), .groups = "drop")
  res$gene = g
  return(res)
})
f3_df = bind_rows(f3_list)

fig3 = ggplot(f3_df, aes(x = sound, y = pct, fill = group)) +
  geom_col(position = "dodge") +
  facet_wrap(~ gene, scales = "free_y")

print(fig3)



sum(rownames(tunga) %in% final_genes)
cells = subset(tunga, features = final_genes)

Idents(cells) = "sound"
core_data = FindMarkers(cells, ident.1 = "Whine", ident.2 = "Whine-chuck")