# ----------------------------------------------------------------------------------------------
# Bisulfite Sequencing Analysis Workflow (GSE202402)
# Dataset: Xmal-RRBS breast carcinoma pre- and post-chemotherapy
# Objective: Identify Differentially Methylated Loci (DMLs) and Regions (DMRs)
# Input: dss_input.csv (processed from Python)
# ----------------------------------------------------------------------------------------------
# B. Bisulfite Sequencing TF Enrichment Analysis
# ---------------------------------------------------------------------------------------------------------------------
# - **Objective:** Identify TF motifs enriched in Differentially methylated Loci (DMLs)
# - **Challenges Faced:**
#     - High memory usage when using matchMotifs() with large PWM lists
#     - Universalmotif conversion failed due to missing values
#     - Batch processing and background correction resulted in memory usage errors
# - **Optimized Approach:**
#   - Reduce background regions (30%)
#     - Limit to JASPAR CORE TF motifs
#     - Batch-wise motif scanning to prevent memory overflow
#     - Filter motifs using Raw Score > 15 (instead of p-values)
# ------------------------------------------------------------------------------------------------------------------------  

# Set working directory
setwd("D:/My projects_Oct 2023/Bisulfite_seq/Bisulfite Sequencing")

# --------------------------------------------------------------------------------------------------------------------
# Step 1: Install and Load Required pacakges
# --------------------------------------------------------------------------------------------------------------------

# Install necessary Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("JASPAR2022", "TFBSTools", "BSgenome.Hsapiens.UCSC.hg38", 
                       "GenomicRanges", "motifmatchr", "Biostrings", "regioneR",
                       "ggplot2", "Matrix", "GenomeInfoDb"), ask = FALSE)

# Load Required Libraries
library(JASPAR2022)    # TF motifs database
library(TFBSTools)     # Motif conversion
library(BSgenome.Hsapiens.UCSC.hg38)  # Human genome for motif scanning
library(GenomicRanges) # Genomic coordinate operations
library(motifmatchr)   # TF motif scanning
library(Biostrings)    # Sequence manipulation
library(regioneR)      # Background region generation
library(ggplot2)       # Visualization
library(Matrix)        # Sparse matrix operations
library(GenomeInfoDb)  # Chromosome formatting
library(BSgenome.Hsapiens.UCSC.hg38.masked)  # For generating background regions
# -----------------------------------------------------------------------------------------------------------------------------
# Step 2: Load and Preprocess DML data
# ------------------------------------------------------------------------------------------------------------------------------
# Load annotated DMLs
annotated_dml = read.csv("annotated_DMLs.csv")

# Convert DMLs to GRanges format
gr_dml = GRanges(seqnames = as.character(annotated_dml$seqnames),
                 ranges = IRanges (start = as.numeric(annotated_dml$start),
                                   end = as.numeric(annotated_dml$start) + 1)) # DMLs are single CpG sites

# Filter out unplaced/alternative chromosomes
seqlevelsStyle(gr_dml) = "UCSC" # Ensure chromosome names match UCSC format
gr_dml = keepStandardChromosomes(gr_dml, pruning.mode = "coarse")

cat("DML regions loaded:", length(gr_dml), "\n")
# DML regions loaded: 202

#-----------------------------------------------------------------------------------------------------------------------------------
# Step 3: Generate Background Regions (Reduced by 30%)
# ----------------------------------------------------------------------------------------------------------------------------------

set.seed(123) # Ensure reproducibility

# Reduce background size to 30% of 'gr_dml'
background_size = round(length(gr_dml) * 0.3)

# Generate background regions (shuffle genomic coordinates)
background_regions_dml = sample(
  randomizeRegions(gr_dml, genome = BSgenome.Hsapiens.UCSC.hg38.masked,
                   per.chromosome = TRUE,
                   allow.overlaps = TRUE),
  size = background_size)

cat("Backgroun regions generated: ", length(background_regions_dml), "\n")
# Backgroun regions generated:  61

# -------------------------------------------------------------------------------------------------------------
# Step 4: Load JASPAR CORE TF Motifs (High-Confidence Motifs)
# -----------------------------------------------------------------------------------------------------

opts <- list(species = 9606, collection = "CORE", all_versions = FALSE)
pfm_list_core <- getMatrixSet(JASPAR2022, opts)  # Retrieve JASPAR CORE motifs

# Convert PFMatrix to PWMatrix, ensuring valid motifs
pwm_list_core <- lapply(pfm_list_core, function(pfm) {
  if (inherits(pfm, "PFMatrix")) {
    pwm <- tryCatch(toPWM(pfm), error = function(e) NULL)  # Handle errors gracefully
    if (!is.null(pwm) && inherits(pwm, "PWMatrix")) return(pwm)
  }
  return(NULL)  # Skip invalid entries
})

# Remove NULL values from PWM list
pwm_list_core <- Filter(Negate(is.null), pwm_list_core)

# Ensure at least one valid PWM matrix
if (length(pwm_list_core) == 0) stop("No valid PWM matrices found!")

# Convert to PWMatrixList format
pwm_list_core <- do.call(PWMatrixList, pwm_list_core)

cat("Successfully loaded and converted", length(pwm_list_core), "TF motifs.\n")

# ----------------------------------------------------------------------------------------------------------
# Step 5: Extract Sequences from the Genome
# --------------------------------------------------------------------------------------------------------
dml_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_dml)
bg_sequences_dml <- getSeq(BSgenome.Hsapiens.UCSC.hg38, background_regions_dml)

cat("Extracted sequences for DMLs and background regions.\n")

# -------------------------------------------------------------------------------------------------------------
# Step 6: Scan for TF motifs in DMLs (batch processing)
# --------------------------------------------------------------------------------------------------------------
batch_size <- 20  # Adjust batch size to reduce memory usage
pwm_batches <- split(pwm_list_core, ceiling(seq_along(pwm_list_core) / batch_size))

# Function to scan motifs in batches
scan_motifs_in_batches <- function(pwm_batches, target_regions) {
  motif_hits_list <- lapply(pwm_batches, function(batch) {
    tryCatch(
      matchMotifs(batch, target_regions, genome = "hg38", out = "scores"),
      error = function(e) {
        cat("Error in batch processing:", e$message, "\n")
        return(NULL)
      }
    )
  })
  # Remove NULL entries
  motif_hits_list <- Filter(Negate(is.null), motif_hits_list)
  return(motif_hits_list)
}

# Scan for TF motifs in **DMLs only** (NO Background Correction)
motif_hits_dml_list <- scan_motifs_in_batches(pwm_batches, gr_dml)

cat("Successfully scanned for TF motifs in DMLs.\n")

# --------------------------------------------------------------------------------------------------------------------
# Step 7: Extract Significant TF Motifs (Raw Score > 15)
# ---------------------------------------------------------------------------------------------------------------------
motif_scores_dml <- do.call(cbind, lapply(motif_hits_dml_list, motifScores))

# Identify motifs exceeding a **raw score > 15**
significant_motifs_1 <- which(motif_scores_dml > 15, arr.ind = TRUE)

# Extract TF names from JASPAR
motif_names <- colnames(motif_scores_dml)
significant_tf_motifs_1 <- unique(motif_names[significant_motifs_1[, 2]])

cat("\nIdentified", length(significant_tf_motifs_1), "Significant TF Motifs (Raw Score > 15):\n")
print(significant_tf_motifs_1)

# ------------------------------------------------------------------------------------------------------------------------
# Step 8: Visualize Top TF motifs
# -----------------------------------------------------------------------------------------------------------------------

# Count occurrences of each TF
tf_counts <- table(significant_tf_motifs_1)

# Convert to dataframe for visualization
tf_df <- data.frame(TF = names(tf_counts), Count = as.numeric(tf_counts))

# Sort by count
tf_df <- tf_df[order(-tf_df$Count), ]

# Ensure at least 15 TFs for plotting
num_tfs <- min(nrow(tf_df), 15)

# Generate Bar Plot for Top TFs
ggplot(tf_df[1:num_tfs, ], aes(x = reorder(TF, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Transcription Factors Enriched in DMLs (Raw Score > 15)",
       x = "Transcription Factor",
       y = "Enrichment Count")

#--------------------------------END OF TF ENRICHMENT ANALYSIS ON DML -----------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------
# B. Bisulfite Sequencing TF Enrichment Analysis
# ---------------------------------------------------------------------------------------------------------------------
# **Objective:** Identify TF motifs enriched in Differentially methylated Regions (DMRs)
#   - JASPAR2022: To retrieve known TF binding motifs.
#   - MotifMatchr: For scanning DMRs for TF motifs.
#   - Background Correction: Using randomized genomic regions to account for non-specific motif enrichment.
#   - Multiple Testing Correction: Applying Benjamini-Hochberg (BH) correction.

# ------------------------------------------------------------------------------------------------------------------------  
# Step 1: Install & Load Required Pacakges
# --------------------------------------------------------------------------------------------------------------------------
# Install required Bioconductor packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("JASPAR2022", "TFBSTools", "BSgenome.Hsapiens.UCSC.hg38", 
                       "motifmatchr", "regioneR", "ggplot2", "Matrix"), ask = FALSE)

# Load necessary libraries
library(JASPAR2022)     # TF motif database
library(TFBSTools)       # Motif processing
library(BSgenome.Hsapiens.UCSC.hg38)  # Human genome
library(motifmatchr)     # TF motif scanning
library(GenomicRanges)   # Genomic coordinate operations
library(regioneR)        # Background region generation
library(Matrix)          # Sparse matrix operations
library(ggplot2)         # Visualization

#----------------------------------------------------------------------------------------------------------------------------
# Step 2: Load JASPAR TF Motifs
# ----------------------------------------------------------------------------------------------------------------------

# Retrieve TF motifs for Homo sapiens (human)
opts <- list(species = 9606, all_versions = FALSE)
pfm_list <- getMatrixSet(JASPAR2022, opts)

# Convert PFMatrix to PWMatrix for motif scanning
pwm_list <- lapply(pfm_list, function(pfm) {
  if (inherits(pfm, "PFMatrix")) {
    pwm <- tryCatch(toPWM(pfm), error = function(e) NULL)  
    if (!is.null(pwm) && inherits(pwm, "PWMatrix")) return(pwm)  
  }
  return(NULL)  # Skip invalid entries
})

# Remove NULL values from PWM list
pwm_list <- Filter(Negate(is.null), pwm_list)

# Ensure valid motif list
if (length(pwm_list) == 0) stop("No valid PWM motifs found!")

# Convert list to PWMatrixList format
pwm_list <- do.call(PWMatrixList, pwm_list)
cat("Successfully loaded", length(pwm_list), "TF motifs.\n")

# ---------------------------------------------------------------------------------------------------------------------------------
# Step 3: Prepare DMR data
# --------------------------------------------------------------------------------------------------------------------------------------
# Load DMR results (Ensure 'dmr_results.csv' is pre-generated)
dmr_results <- read.csv("significant_DMRs.csv")

# Convert DMRs to GRanges format
gr_dmr <- GRanges(seqnames = as.character(dmr_results$chr),
                  ranges = IRanges(start = as.numeric(dmr_results$start), 
                                   end = as.numeric(dmr_results$end)))

cat("Loaded", length(gr_dmr), "DMRs for analysis.\n")

# ----------------------------------------------------------------------------------------------------------------
# Generate Background regions for Enrichment Analysis
# -----------------------------------------------------------------------------------------------------------------
# Ensure masked genome is installed for proper background correction
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38.masked", force = TRUE, ask = FALSE)
library(BSgenome.Hsapiens.UCSC.hg38.masked)

# Set random seed for reproducibility
set.seed(123)

# Generate **randomized background regions** matching DMR size distribution
background_regions <- randomizeRegions(gr_dmr, 
                                       genome = BSgenome.Hsapiens.UCSC.hg38.masked, 
                                       per.chromosome = TRUE)

cat("Generated", length(background_regions), "background regions.\n")

# ------------------------------------------------------------------------------------------------------------------------
# Step 5: Extract DNA sequences
# ----------------------------------------------------------------------------------------------------------------------
# Extract sequences from DMRs & Background Regions
dmr_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_dmr)
bg_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, background_regions)

cat("Extracted sequences for DMRs and background.\n")

#---------------------------------------------------------------------------------------------------------------------------------
# Scan for TF motifs in DMR and background
# ----------------------------------------------------------------------------------------------------------
# Scan TF motifs in DMRs
motif_hits_dmr <- matchMotifs(pwm_list, gr_dmr, genome = "hg38", out = "scores")

# Scan TF motifs in background regions
motif_hits_bg <- matchMotifs(pwm_list, background_regions, genome = "hg38", out = "scores")

cat("TF motif scanning completed.\n")

# ------------------------------------------------------------------------------------------------------------------------------
# Step 7: Compute Enrichment Scores
# --------------------------------------------------------------------------------------------------------
# Compute motif enrichment (Fold Enrichment)
motif_enrichment <- rowMeans(motifScores(motif_hits_dmr) + 1) / rowMeans(motifScores(motif_hits_bg) + 1)

# Extract motif scores
motif_scores <- motifScores(motif_hits_dmr)

# Identify motifs with raw score > 15 (High-confidence motifs)
significant_motifs_1 <- rownames(motif_scores)[apply(motif_scores, 1, max) > 15]

cat("Identified", length(significant_motifs_1), "motifs with raw score > 15.\n")

# -------------------------------------------------------------------------------------------------------------------------------
# Step 8: Statistical Significance Testing
# ---------------------------------------------------------------------------------------------------------
# Convert log-scores to approximate p-values
p_values <- 10^(-motif_scores)

# Adjust p-values using Benjamini-Hochberg (BH) correction
adjusted_p <- p.adjust(apply(p_values, 1, min), method = "BH")

# Identify motifs based on **adjusted p-value < 0.05** and **log2FC > 1**
significant_motifs_2 <- rownames(motif_scores)[(motif_enrichment > 1) & (adjusted_p < 0.05)]

cat("Identified", length(significant_motifs_2), "significant motifs (FDR < 0.05).\n")

# -------------------------------------------------------------------------------------------------------
# Step 9: Extract and save significant TF motifs
# ---------------------------------------------------------------------------------------------------------------
# Extract TF names from JASPAR
motif_names <- colnames(motif_scores)

# Get significant motifs
significant_tf_motifs_1 <- unique(motif_names[significant_motifs_1])
significant_tf_motifs_2 <- unique(motif_names[significant_motifs_2])

# Save results
write.csv(significant_tf_motifs_1, "TF_enrichment_rawscore.csv", row.names = FALSE)
write.csv(significant_tf_motifs_2, "TF_enrichment_FDR.csv", row.names = FALSE)

cat("Saved TF enrichment results.\n")

# --------------------------------------------------------------------------------------------------------------
# Step 10: Visualization: Top Enriched Transcription Factors
# --------------------------------------------------------------------------------------------------------
# Count occurrences of each TF
tf_counts <- table(significant_tf_motifs_1)

# Convert to dataframe for visualization
tf_df <- data.frame(TF = names(tf_counts), Count = as.numeric(tf_counts))

# Sort by count
tf_df <- tf_df[order(-tf_df$Count), ]

# Ensure at least 15 TFs for plotting
num_tfs <- min(nrow(tf_df), 15)

# Generate Bar Plot for Top TFs
ggplot(tf_df[1:num_tfs, ], aes(x = reorder(TF, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Enriched Transcription Factors in DMRs",
       x = "Transcription Factor",
       y = "Enrichment Count")

ggsave("TF_enrichment_plot.png", width = 8, height = 6, dpi = 300)
cat("TF enrichment plot saved.\n")

# -------------------------------END OF TF ENRICHMENT ANALYSIS ON DMR -------------------------------------------------------------------

