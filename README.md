# Tungara Frog Neural snRNA-seq

## About the Project
This project analyzes the brains of Tungara frogs. The research is conducted in the Dr. Hoffman Lab. The goal is to potentially understand how different mating calls ("Whine" vs. "Whine-Chuck") might affect neural activity. Single-nucleus RNA sequencing (snRNA-seq) is utilized to study this at the cellular level, which could possibly help model theoretical neural response circuits. Transcriptomic visualizations were generated, which contributed to a third authorship on a research poster presented at the 2026 Joint Meeting of Ichthyologists and Herpetologists.

## Biological Implications (Preliminary)
The data suggests the following possible neural circuitry model using key cell roles and pathways:
- Location: The response might originate in the Torus Semicircularis (auditory midbrain), potentially bounded by specific developmental genes.
- Trigger: The complex "Whine-Chuck" call could potentially bypass baseline neural blockers.  Senders: "Factory" cells might use early genes to deploy signal ligands.
- Receivers: "Receiver" cells may bind these ligands, which could activate internal memory genes (CAMK1, CREB1, MAPK8) to tentatively unlock the genome.
- Two Possible Pathways: The response might then split into two possible strategies:Endocrine Path: Could be driven by CAPRIN2 to pump hormones and neuropeptides.
- Structural Path: Could be driven by NRTN to rewire cell synapses.  Genetics: The choice between Endocrine or Structural paths might depend heavily on the individual frog's genotype.
- Action: "Motor" cells might act as the final step to deploy axon guidance and potentially trigger movement.  

## Tools and Methods
Because the Tungara frog lacks a fully annotated genome, custom tools had to be built to map the data.
- Genome Translation: A custom BLASTp Reciprocal Best Hit pipeline was written in Bash. This utilizes RefSeq protein files for forward and reverse BLAST searches, isolating 1-to-1 true orthologs to translate the Tungara genome using the Xenopus reference.
- Demultiplexing & QC: Souporcell SNP calling was utilized to attempt to demultiplex the data to individual frogs. Pseudobulk matrices were built to filter out suspected dropout-driven fold changes and possible pseudoreplication artifacts.
- Data Processing: The Seurat package in R was used to clean, normalize, and integrate data via CCA, which could possibly correct for batch effects.
- Dimensionality Reduction: UMAP and PCA were utilized to visualize possible brain cell types.
- Cell Annotation & Validation: Sc-type was used for automated labeling paired with custom BLAST dictionaries. Clustering parameters were compared against an independent pipeline using Adjusted Rand Index (ARI) scores and alluvial diagrams.
- Network Modeling: hdWGCNA was used to build co-expression networks and find potential "hub genes". Broad modules were fractured to guess at specific functional cellular strategies.
- Neural Circuit Mapping: Theoretical circuits were tested by mapping the co-expression of regulatory genes (CAPRIN2, NRTN, MAPK8) via correlation matrices and functional enrichment.  

## Repository Structure
- scripts/blastp_pipeline.sh: The Bash script for cross-species genome translation via two-way BLAST searches.
- scripts/seurat_processing.R: The R code for data cleaning, CCA integration, and UMAP clustering.
- scripts/hdwgcna_network.R: The R code for building gene co-expression networks.
- hedi_frog_sctype_method.R: The script for automated cell annotation and differential expression testing.
- hedi_frog_by_blast.R: The script integrating the BLASTp translation table and setting up hdWGCNA.
- method_comparaison.R: The pipeline validation script comparing cell labels via Adjusted Rand Index scoring.
- biological_implication_justification.R: The script attempting to model the auditory circuit from baseline repression to motor signal.
- data_genotypes.R: The script for demultiplexing data to individual frogs and applying pseudobulk QC.
- hedi_frog_finalanalysis.R: The script fracturing active gene modules to explore Endocrine vs. Structural strategies mapped to genotypes.

- 
