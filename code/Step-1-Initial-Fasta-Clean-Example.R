# R 2025.05.1+513
# 8th June 2026
# Reviewed

# ── Dependencies ──────────────────────────────────────────────────────────────
library(Biostrings)
library(ape)
library(ggtree)
library(ggplot2)


# ── 1. Load and deduplicate raw GISAID sequences ──────────────────────────────

fasta <- readDNAStringSet("gisaid_epiflu_sequence.fasta")

# GISAID downloads often contain exact duplicates — remove any sequences
# where both the header *and* the sequence string are identical
fasta_uniq <- fasta[!duplicated(names(fasta)) & !duplicated(as.character(fasta))]

writeXStringSet(fasta_uniq, "gisaid_epiflu_sequence_2.fasta")


# ── 2. Filter out unwanted sequences ──────────────────────────────────────────

fasta <- readDNAStringSet("gisaid_epiflu_sequence_2.fasta")

# Drop sequences flagged as embargoed or from faecal samples — these are either
# unavailable for sharing or generally unreliable for HA phylogenetics
exclude_terms <- c("feces", "embargo")
pattern       <- paste(exclude_terms, collapse = "|")
keep          <- !grepl(pattern, names(fasta), ignore.case = TRUE)
fasta_filtered <- fasta[keep]


# ── 3. Clean up FASTA headers ─────────────────────────────────────────────────

# Replace characters that break downstream tools (MAFFT, IQ-TREE) with underscores
names(fasta_filtered) <- gsub("[()']", "_", names(fasta_filtered))

# Remove "HA|" segment prefix — inherited from GISAID multisegment exports,
# not needed once we're working with HA only
names(fasta_filtered) <- gsub("HA\\|", "", names(fasta_filtered))

head(names(fasta_filtered))  # quick sanity check before writing out

writeXStringSet(fasta_filtered, "input.fasta")


# ── 4. Multiple sequence alignment ───────────────────────────────────────────

# --auto lets MAFFT choose the appropriate strategy based on dataset size;
# fine for most HA datasets, but swap for --localpair / --maxiterate 1000
# if you want higher accuracy on divergent sequences
system("mafft --auto input.fasta > aligned.fasta")


# ── 5. Quality filtering: remove gappy / ambiguous sequences ─────────────────

aligned <- readDNAStringSet("aligned.fasta")

# Calculate per-sequence proportion of gaps and N's separately, then combine.
# High gap/N content usually indicates poor-quality or short sequences that
# would distort branch lengths in the tree.
gap_freq      <- letterFrequency(aligned, letters = "-", as.prob = TRUE)
n_freq        <- letterFrequency(aligned, letters = "N", as.prob = TRUE)
combined_freq <- gap_freq + n_freq

cat("Gap + N frequency summary:\n")
print(summary(as.numeric(combined_freq)))

# Visual distribution to inform threshold choice
hist(
  combined_freq,
  breaks = 50,
  main   = "Gap + N frequency per sequence",
  xlab   = "Proportion of gaps (-) and ambiguous bases (N)",
  col    = "steelblue"
)

# !! Adjust this threshold based on the histogram above !!
# 0.3 is a reasonable starting point for HA; tighter datasets may warrant 0.1–0.2
threshold     <- 0.3
keep          <- combined_freq < threshold
aligned_clean <- aligned[keep]

writeXStringSet(aligned_clean, "aligned_clean.fasta")


# ── 6. Manual trimming step ───────────────────────────────────────────────────

# Inspect aligned_clean.fasta in, e.g., Geneious or AliView and trim ragged
# ends / UTR regions before proceeding. Save the trimmed file as
# aligned_clean_trim.fasta (used later in the lineage-splitting step).


# ── 7. Re-align after trimming, then build ML tree ───────────────────────────

# Re-align after any manual edits to ensure consistency
system("mafft --auto aligned_clean.fasta > aligned_clean_2.fasta")

# GTR+G is standard for influenza HA; 1000 ultrafast bootstraps for node support.
# -nt AUTO lets IQ-TREE pick thread count — set a fixed integer on shared servers
# to avoid monopolising cores.
system("iqtree2 -redo -s aligned_clean_2.fasta -m GTR+G -bb 1000 -nt AUTO -pre iqtree_out")


# ══════════════════════════════════════════════════════════════════════════════
#  OPTIONAL: Split dataset by lineage (e.g. H3 Avian Eurasian vs Americas)
#  Only needed when your alignment spans multiple distinct lineages.
#  Workflow: open iqtree_out.treefile in FigTree → manually define clade
#  selections → export as iqtree_out_lineages.tree (Newick with clade labels).
# ══════════════════════════════════════════════════════════════════════════════

# ── 8. Extract lineage tip labels from annotated tree ────────────────────────

tree_lin <- read.tree("iqtree_out_lineages.tree")

# Tip label vectors come from named subtrees defined in FigTree.
# The list element names ("treetree_7=", "treetree_8=") are FigTree's internal
# clade identifiers — check these against your exported .tree file if they differ.
eurasian_lineage_tips <- tree_lin[["treetree_7="]][["tip.label"]]
americas_lineage_tips <- tree_lin[["treetree_8="]][["tip.label"]]

# FigTree sometimes wraps tip labels in single quotes — strip them so they
# match the FASTA headers exactly
eurasian_lineage_tips <- gsub("^'|'$", "", eurasian_lineage_tips)
americas_lineage_tips <- gsub("^'|'$", "", americas_lineage_tips)


# ── 9. Subset FASTA by lineage and write out ──────────────────────────────────

# Work from the trimmed alignment so both output files are analysis-ready
fasta <- readDNAStringSet("aligned_clean_trim.fasta")

eurasian_fasta <- fasta[names(fasta) %in% eurasian_lineage_tips]
americas_fasta <- fasta[names(fasta) %in% americas_lineage_tips]

cat("Eurasian sequences:", length(eurasian_fasta), "\n")
cat("Americas sequences:", length(americas_fasta), "\n")

writeXStringSet(eurasian_fasta, "H3_Eu.fasta")
writeXStringSet(americas_fasta, "H3_Am.fasta")

# Each output file can now be taken through steps 4–7 independently to
# produce lineage-specific trees.
