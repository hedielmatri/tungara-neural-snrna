#!/bin/bash

WORKDIR="/stor/scratch/FRI-BigDataBio/FRI_summer_2026/frog_data/Hedi_R"
cd $WORKDIR
THREADS=32

echo "Tungara-Xenopus 1-to-1 Ortholog Pipeline"

# DOWNLOAD THE REFSEQ FILES
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/894/005/GCF_040894005.1_aEngPut4.maternal/GCF_040894005.1_aEngPut4.maternal_protein.faa.gz
gunzip -f GCF_040894005.1_aEngPut4.maternal_protein.faa.gz

# Xenopus tropicalis (UCB_Xtro_10.0)
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.4_UCB_Xtro_10.0/GCF_000004195.4_UCB_Xtro_10.0_protein.faa.gz
gunzip -f GCF_000004195.4_UCB_Xtro_10.0_protein.faa.gz

TUNGARA_PROT="GCF_040894005.1_aEngPut4.maternal_protein.faa"
REF_PROT="GCF_000004195.4_UCB_Xtro_10.0_protein.faa"
# TWO-WAY BLAST

makeblastdb -in $TUNGARA_PROT -dbtype prot -out tungara_refseq_db
makeblastdb -in $REF_PROT -dbtype prot -out xenopus_refseq_db

echo "      Running Forward BLAST (Tungara -> Xenopus)..."
blastp -query $TUNGARA_PROT -db xenopus_refseq_db -out tungara_vs_ref.blastp.tsv \
  -evalue 1e-5 -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" -num_threads $THREADS

echo "      Running Reverse BLAST (Xenopus -> Tungara)..."
blastp -query $REF_PROT -db tungara_refseq_db -out ref_vs_tungara.blastp.tsv \
  -evalue 1e-5 -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" -num_threads $THREADS

# RECIPROCAL BEST HITS

# Tungara -> Xenopus
awk '{
  q = $1; bits = $8 + 0
  if (!(q in best_bits) || bits > best_bits[q]) {
    best_bits[q] = bits; best_line[q] = $0
  }
} END { for (q in best_line) print best_line[q] }' tungara_vs_ref.blastp.tsv > tungara_vs_ref.best.tsv

# Xenopus -> Tungara
awk '{
  q = $1; bits = $8 + 0
  if (!(q in best_bits) || bits > best_bits[q]) {
    best_bits[q] = bits; best_line[q] = $0
  }
} END { for (q in best_line) print best_line[q] }' ref_vs_tungara.blastp.tsv > ref_vs_tungara.best.tsv

# Intersections
echo "xenopus_protein,tungara_protein" > xenopus_to_tungara_RBH.csv
awk '
NR==FNR { tun_to_ref[$1] = $2; next }
{ ref_to_tun[$1] = $2 }
END {
  count = 0
  for (tun in tun_to_ref) {
    ref = tun_to_ref[tun]
    if ((ref in ref_to_tun) && ref_to_tun[ref] == tun) {
      print ref "," tun
      count++
    }
  }
  print "Found " count " 1-to-1 true orthologs" > "/dev/stderr"
}' tungara_vs_ref.best.tsv ref_vs_tungara.best.tsv >> xenopus_to_tungara_RBH.csv

# 4. EXTRACT DICTIONARIES
# --------------------------------------------------

wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/654/675/GCF_017654675.1_Xenopus_laevis_v10.1/GCF_017654675.1_Xenopus_laevis_v10.1_feature_table.txt.gz
echo "protein_id,gene_symbol" > xenopus_id_to_symbol.csv
zcat GCF_017654675.1_Xenopus_laevis_v10.1_feature_table.txt.gz | awk -F'\t' '
NR==1 { for(i=1; i<=NF; i++) { if($i=="product_accession") p=i; if($i=="symbol") s=i; if($i=="GeneID") g=i } }
NR>1 && $p != "" { sym = $s != "" ? $s : "LOC" $g; if(sym != "") print $p "," sym }
' >> xenopus_id_to_symbol.csv

wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/040/894/005/GCF_040894005.1_aEngPut4.maternal/GCF_040894005.1_aEngPut4.maternal_genomic.gff.gz
echo "protein_id,gene_symbol" > tungara_id_to_symbol_FIXED.csv
zcat GCF_040894005.1_aEngPut4.maternal_genomic.gff.gz | awk -F'\t' '$3=="CDS" {print $9}' | awk '{
  prot=""; gene="";
  if (match($0, /[XN]P_[0-9]+\.[0-9]+/)) { prot = substr($0, RSTART, RLENGTH); }
  if (match($0, /gene=[^;]+/)) { gene = substr($0, RSTART+5, RLENGTH-5); }
  if (prot != "" && gene != "") print prot "," toupper(gene);
}' | sort -u >> tungara_id_to_symbol_FIXED.csv


# tropicalis
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.4_UCB_Xtro_10.0/GCF_000004195.4_UCB_Xtro_10.0_genomic.gff.gz

echo "protein_id,gene_symbol" > xenopus_id_to_symbol_tropicalis.csv

zcat GCF_000004195.4_UCB_Xtro_10.0_genomic.gff.gz | awk -F'\t' '$3=="CDS" {print $9}' | awk '{
  prot=""; gene="";
  if (match($0, /[XN]P_[0-9]+\.[0-9]+/)) { prot = substr($0, RSTART, RLENGTH); }
  if (match($0, /gene=[^;]+/)) { gene = substr($0, RSTART+5, RLENGTH-5); }
  if (prot != "" && gene != "") print prot "," toupper(gene);
}' | sort -u >> xenopus_id_to_symbol_tropicalis.csv
