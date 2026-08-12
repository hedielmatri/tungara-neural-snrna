setwd("/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/Hedi_R")

library(Seurat)
library(dplyr)
library(ggplot2)

options(future.globals.maxSize = 64000 * 1024^2)
Sys.setenv(R_MAX_VSIZE = "64Gb")

whole_brain = readRDS("neural_wgcna_obj.rds")

whole_brain = subset(whole_brain, orig.ident != "Green noise")
whole_brain$sound_type = ifelse(whole_brain$orig.ident == "Whine", "Baseline", "Complex")

# NRTN CAPRIN2: "NRTN acts as the master plasticity output while CAPRIN2 binds and stabilizes its newly minted transcripts"
if (!file.exists("core_data.rds")) {
  Idents(whole_brain) = "sound_type"
  core_data = FindMarkers(whole_brain, ident.1 = "Complex", ident.2 = "Baseline")
  saveRDS(core_data, "core_data.rds")
}
core_data = readRDS("core_data.rds")
core_plot = VlnPlot(whole_brain, features = c("NRTN", "CAPRIN2"))
print(core_plot)


cell_labels = whole_brain[["refined_annotation"]][, 1]
whole_brain$is_epicenter = grepl("sst|tac1", cell_labels)
whole_brain$brain_region = ifelse(whole_brain$is_epicenter, "Epicenter", "Other")
Idents(whole_brain) = "brain_region"

# PAX2 PAX3 PAX7 GATA2 GATA3: "These developmental regulators define the Torus Semicircularis auditory midbrain boundaries"
if (!file.exists("anatomy_data.rds")) {
  anatomy_data = FindMarkers(whole_brain, ident.1 = "Epicenter", ident.2 = "Other")
  saveRDS(anatomy_data, "anatomy_data.rds")
}
anatomy_data = readRDS("anatomy_data.rds")
anatomy_plot = DotPlot(whole_brain, features = c("PAX2", "PAX3", "PAX7", "GATA2", "GATA3"))
# Check the dots to see the local brain area
print(anatomy_plot)


epicenter = subset(whole_brain, idents = "Epicenter")
Idents(epicenter) = "sound_type"

# BOC LOC140111485: "These genes enforce a baseline blockade to repress axon guidance and prevent motor projection"
if (!file.exists("baseline_data.rds")) {
  baseline_data = FindMarkers(epicenter, ident.1 = "Baseline", ident.2 = "Complex")
  saveRDS(baseline_data, "baseline_data.rds")
}
baseline_data = readRDS("baseline_data.rds")
baseline_plot = VlnPlot(epicenter, features = c("BOC", "LOC140111485"))
print(baseline_plot)


complex_sound = subset(epicenter, idents = "Complex")
complex_labels = complex_sound[["refined_annotation"]][, 1]

complex_sound$is_factory = grepl("sst", complex_labels)
complex_sound$cell_role = ifelse(complex_sound$is_factory, "Factory", "Receiver")
Idents(complex_sound) = "cell_role"

# EGR3 EGR4 PLAT: "EGR3 and EGR4 early genes that upregulate to flood the circuit with PLAT synaptic cleavers"
if (!file.exists("factory_data.rds")) {
  factory_data = FindMarkers(complex_sound, ident.1 = "Factory", ident.2 = "Receiver")
  saveRDS(factory_data, "factory_data.rds")
}
factory_data = readRDS("factory_data.rds")
factory_plot = DotPlot(complex_sound, features = c("EGR3", "EGR4", "PLAT"))
print(factory_plot)


# FLRT1: "Ligand deployed by tripwire cells to talk  with receiver cells"
ligand_plot = VlnPlot(complex_sound, features = "FLRT1")
# Spot the physical signal between the cells
print(ligand_plot)


# SLC6A5: "An inhibitory glycine transporter used by receiver cells to filter out background noise"
filter_plot = VlnPlot(complex_sound, features = "SLC6A5")
print(filter_plot)


# test ADGRA3 mathematically independently across sounds
all_receivers = subset(epicenter, grepl("tac1", refined_annotation))
Idents(all_receivers) = "sound_type"

# ADGRA3: "The adhesion receptor that binds the physical ligand to trigger a localized calcium influx"
if (!file.exists("receptor_data.rds")) {
  receptor_data = FindMarkers(all_receivers, ident.1 = "Complex", ident.2 = "Baseline")
  saveRDS(receptor_data, "receptor_data.rds")
}
receptor_data = readRDS("receptor_data.rds")
receptor_plot = VlnPlot(all_receivers, features = "ADGRA3")
print(receptor_plot)


# CAMK1 CREB1 MAPK8 ZNF292: "CAMK1 triggers CREB1 and MAPK8 internal relays which unlock the genome via ZNF292 zinc finger transcription factors"
wire_plot = DotPlot(all_receivers, features = c("CAMK1", "CREB1", "MAPK8", "ZNF292"))
print(wire_plot)


# SLC22A2 DDAH1 HIP1R WIF1: "SLC22A2 shuts off to trap dopamine while HIP1R swallows synaptic anchors and DDAH1 blasts nitric oxide and WIF1 blocks structural growth"
lockdown_plot = VlnPlot(all_receivers, features = c("SLC22A2", "DDAH1", "HIP1R", "WIF1"))
print(lockdown_plot)


complex_sound$is_motor = grepl("Excitatory", complex_labels)
complex_sound$motor_role = ifelse(complex_sound$is_motor, "Motor", "Filter")
Idents(complex_sound) = "motor_role"

# SYT7 SYT9 EPHA8 PLXNA1: "SYT7 and SYT9 act as exocytosis sensors that deploy EPHA8 and PLXNA1 axon guidance pathways to punch out the final motor signal"
if (!file.exists("motor_data.rds")) {
  motor_data = FindMarkers(complex_sound, ident.1 = "Motor", ident.2 = "Filter")
  saveRDS(motor_data, "motor_data.rds")
}
motor_data = readRDS("motor_data.rds")
exit_plot = DotPlot(complex_sound, features = c("SYT7", "SYT9", "EPHA8", "PLXNA1"))
print(exit_plot)

