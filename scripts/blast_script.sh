#!/bin/bash

WORKDIR="/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/Hedi_R"
cd $WORKDIR

# Define your input files (Ensure these match your downloaded NCBI files)
TUNGARA_PROT="GCA_019512145.1_UCB_Epus_1.0_protein.faa"
REF_PROT="GCF_017654675.1_Xenopus_laevis_v10.1_protein.faa"
THREADS=32

# 1. DATABASES
# --------------------------------------------------

makeblastdb \
  -in $TUNGARA_PROT \
  -dbtype prot \
  -out tungara_db

makeblastdb \
  -in $REF_PROT \
  -dbtype prot \
  -out reference_db

# 2. RUN BLAST IN BOTH DIRECTIONS 
# --------------------------------------------------

blastp \
  -query $TUNGARA_PROT \
  -db reference_db \
  -out tungara_vs_ref.blastp.tsv \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  -num_threads $THREADS

blastp \
  -query $REF_PROT \
  -db tungara_db \
  -out ref_vs_tungara.blastp.tsv \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" \
  -num_threads $THREADS

# 3. FIND THE BEST HITS FOR BOTH DIRECTIONS
# --------------------------------------------------

awk '{
  q = $1; bits = $8 + 0
  if (!(q in best_bits) || bits > best_bits[q]) {
    best_bits[q] = bits; best_line[q] = $0
  }
} END { for (q in best_line) print best_line[q] }' tungara_vs_ref.blastp.tsv > tungara_vs_ref.best.tsv

awk '{
  q = $1; bits = $8 + 0
  if (!(q in best_bits) || bits > best_bits[q]) {
    best_bits[q] = bits; best_line[q] = $0
  }
} END { for (q in best_line) print best_line[q] }' ref_vs_tungara.blastp.tsv > ref_vs_tungara.best.tsv

# 4. FIND RECIPROCAL BEST HITS 
# --------------------------------------------------

echo "xenopus_gene,tungara_gene" > xenopus_gene_to_tungara_gene.csv

# Find the 1-to-1 matches and append them as comma-separated values
awk '
NR==FNR { tun_to_ref[$1] = $2; next }
{ ref_to_tun[$1] = $2 }
END {
  for (tun in tun_to_ref) {
    ref = tun_to_ref[tun]
    if ((ref in ref_to_tun) && ref_to_tun[ref] == tun) {
      print ref "," tun
    }
  }
}' tungara_vs_ref.best.tsv ref_vs_tungara.best.tsv >> xenopus_gene_to_tungara_gene.csv

