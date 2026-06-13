# Tungara Frog Neural snRNA-seq Profiling

## About the Project
This project analyzes the brains of Tungara frogs. I am doing this research in the Dr. Hoffman Lab. The goal is to understand how different mating calls ("Whine" vs. "Whine-Chuck") affect the frog's neural activity. I am using single-nucleus RNA sequencing (snRNA-seq) to study this at the cellular level.

## Current Status: Active Research
This project is currently ongoing. I have built the main data processing pipelines. I have successfully translated the genome, clustered the cells, and built the co-expression networks. I am currently running differential expression tests to find specific genes linked to auditory processing. Final biological conclusions will be added when the research is complete.

## Tools and Methods
Because the Tungara frog does not have a fully annotated genome, I had to build custom tools to map the data.
* **Genome Translation:** I wrote a custom BLASTp Reciprocal Best Hit pipeline in Bash. This translated the Tungara genome using the Xenopus (African clawed frog) reference genome.
* **Data Processing:** I used the Seurat package in R to clean, normalize, and integrate the data using CCA.
* **Dimensionality Reduction:** I used UMAP and PCA to visualize the different types of brain cells.
* **Cell Annotation:** I used `sc-type` for automated cell type labeling (e.g., Excitatory Neurons, GABAergic Neurons).
* **Network Modeling:** I used `hdWGCNA` to build co-expression networks and find "hub genes" related to the mating calls.

## Repository Structure
* `scripts/blastp_pipeline.sh`: My Bash script for cross-species genome translation.
* `scripts/seurat_processing.R`: The R code for data cleaning, CCA integration, and UMAP clustering.
* `scripts/hdwgcna_network.R`: The R code for building the gene co-expression networks.

