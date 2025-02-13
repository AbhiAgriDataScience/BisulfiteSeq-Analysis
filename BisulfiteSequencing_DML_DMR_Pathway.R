# ----------------------------------------------------------------------------------------------
# Bisulfite Sequencing Analysis Workflow (GSE202402)
# Dataset: Xmal-RRBS breast carcinoma pre- and post-chemotherapy
# Objective: Identify Differentially Methylated Loci (DMLs) and Regions (DMRs)
# Input: dss_input.csv (processed from Python)
# ----------------------------------------------------------------------------------------------
# A. DML Analysis with Pathway Enrichment
# -------------------------------------------------------------------------------------------------
# Step 1: Install and load required packages
# ------------------------------------------------------------------------------------------------
# Install necessary Bioconductor packages if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("DSS", "bsseq", "ggplot2", "reshape2", "pheatmap", "coin", 
                       "GenomicRanges", "biomaRt", "ChIPseeker", "clusterProfiler",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"), ask = FALSE)

# Load Required Libraries
library(DSS)           # DML & DMR analysis
library(bsseq)         # Methylation analysis
library(ggplot2)       # Visualization
library(reshape2)      # Data reshaping
library(pheatmap)      # Heatmaps
library(coin)          # Statistical tests
library(GenomicRanges) # Genomic operations
library(biomaRt)       # Gene annotation
library(ChIPseeker)    # Peak annotation
library(clusterProfiler) # Pathway enrichment
library(org.Hs.eg.db)  # Gene ontology mapping
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # hg38 genome annotation

# Set working directory
setwd("D:/My projects_Oct 2023/Bisulfite_seq/Bisulfite Sequencing")

# ------------------------------------------------------------------------------------------------------
# Step 2: Load and preprocess methylation data
# ------------------------------------------------------------------------------------------------------

# Load DSS input file (methylation counts for CpG sites
dss_data = read.csv("dss_input.csv")

# Ensure numeric values for methylation counts
dss_data$Total_C = as.numeric(dss_data$Total_C)
dss_data$Methylated_C = as.numeric(dss_data$Methylated_C)

# Convert methylation data to BSseq format
# - BSseq object is the required input format for DSS (Dispersion Shrinkage for Sequencing) data
# - DSS is used for differntial methylation analysis, particularly in bisulfite sequencing datasets

# Split data by sample ID
# - # Each element in list is a df containing cpG site methylation counts for a single Sample_ID
sample_list = split(dss_data, dss_data$Sample_ID) 

# Transform each sample into required BSseq format
bsseq_input = lapply(sample_list, function(df){
  data.frame(chr = df$Chromosome, pos = df$Position, N = df$Total_C, X = df$Methylated_C) # for each sample, we reformat the data into a simplified structure
})
# The new structure has 4 required columns: chr, pos, N (total cytosines observed/ coverage), X (methylated cytosines observed)

# Extract sample names
sample_names = unique(dss_data$Sample_ID)

# Create BSseq object
# - A BSseq object is a special format used for bisulfite sequencing data in R.
# - It allows for efficient storage, easy access and differential methylation analysis with DSS
# - makesBSseqData() converts the list of methylation df (bsseq_input) into BSseq object, required for DSS analysis
# - sampleNames argument links each dataset to the corresponding Sample_ID
BSobj = makeBSseqData(bsseq_input, sampleNames = sample_names)

# ------------------------------------------------------------------------------------------------------------------------------
# Step 3: Assign Sample metadata and filter low-coverage sites
# ------------------------------------------------------------------------------------------------------------------------

#Define sample conditions (Pre-and post-chemotherapy)
design_matrix = data.frame(
  Sample_ID = sampleNames(BSobj),
  Condition = c("Pre-chemo", "Post_chemo", "Pre-chemo", "Pre-chemo", "Post_chemo", "Post_chemo")
)
rownames(design_matrix) = sampleNames(BSobj)
pData(BSobj) = design_matrix

# Filter out CpG sites with low coverage (CpG sites with atleast 3 reads)
BSobj_filtered = BSobj[rowSums(getCoverage(BSobj)) > 3, ]

# -------------------------------------------------------------------------------------------------------------------------
# Step 4: Perform Differential Methylation Loci Analysis (DML)
# -------------------------------------------------------------------------------------------------------------------------
# 4.1 Perform DML analysis using DSS
# Fit dispersion model
dml_fit = DMLfit.multiFactor(BSobj_filtered, design = pData(BSobj), formula = ~Condition)

# Perform DML testing
dml_test = DMLtest.multiFactor(dml_fit, coef = "ConditionPre-chemo")

# Extract significant DMLs (FDR < 0.05 & abs(stat) > 2)
significant_dml = subset(dml_test, fdrs < 0.05 & abs(stat) > 2)
significant_dml$site = paste(significant_dml$chr, ":", significant_dml$pos)

# Save results 
write.csv(significant_dml, "significant_DMLs.csv", row.names = FALSE)
cat("Significant DMLs saved (",nrow(significant_dml),"loci)\n")

# Significant DMLs saved ( 205 loci)

# 4.2: Wilcoxon Rank-Sum test for validation
# Extract methylation values
meth_values <- getMeth(BSobj_filtered, type = "raw")

# Assign CpG site identifiers
rownames(meth_values) <- paste0(seqnames(rowRanges(BSobj_filtered)), ":", start(rowRanges(BSobj_filtered)))

# Filter High-Variance CpG sites (Variance > 0.01)
site_variances <- apply(meth_values, 1, function(x) var(x, na.rm = TRUE))
valid_sites <- rownames(meth_values)[site_variances > 0.01]
valid_sites <- na.omit(trimws(valid_sites))

# Subset methylation values for common sites
common_sites <- intersect(rownames(meth_values), valid_sites)
filtered_meth_values <- meth_values[common_sites, , drop = FALSE]

# Wilcoxon Rank-Sum Test
group1 <- which(pData(BSobj_filtered)$Condition == "Pre_chemo")
group2 <- which(pData(BSobj_filtered)$Condition == "Post_chemo")

wilcox_results <- apply(filtered_meth_values, 1, function(x) {
  if (length(unique(x[group1])) > 1 & length(unique(x[group2])) > 1) {
    return(wilcox.test(x[group1], x[group2], exact = FALSE)$p.value)
  } else {
    return(NA)
  }
})

# Adjust p-values
adjusted_pvalues <- p.adjust(wilcox_results, method = "BH")

# Extract significant CpG sites
wilcox_df <- data.frame(
  Site = names(wilcox_results),
  P_Value = wilcox_results,
  Adjusted_P_Value = adjusted_pvalues
)
significant_wilcox <- subset(wilcox_df, Adjusted_P_Value < 0.05)

# No signficant sites were observed
# Save results
write.csv(significant_wilcox, "wilcoxon_significant_sites.csv", row.names = FALSE)

# ------------------------------------------------------------------------------------------------------------------------
# Step 5: Functional annotation of DMLs
# ------------------------------------------------------------------------------------------------------------------------

# Convert Significant DMLs to GRanges format
gr_dml = GRanges(seqnames = significant_dml$chr,
                 ranges = IRanges(start = significant_dml$pos,
                                  end = significant_dml$pos)
                 )
# - gr_dml is now a GRange object containing genomic coordinates of DMLs.

library(GenomeInfoDb)
# Ensure chromosome names match UCSC format
seqlevelsStyle(gr_dml) = "UCSC"

# Keep only standard chromosomes (chr 1-22, X, Y, M)
gr_dml = keepStandardChromosomes(gr_dml, pruning.mode = "coarse")

# Annotate DML sites to genomic features (genes, pronmoters, exons, etc.)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

peak_annotation = annotatePeak(gr_dml, TxDb = txdb, tssRegion = c(-2000, 2000), verbose = FALSE) # specifies a window of 2 kb upstream and downstream of TSS to annotate promoter regions

annotated_dml = as.data.frame(peak_annotation)

# Save annotated results
write.csv(peak_annotation, "annotated_DMLs.csv", row.names = FALSE)
cat("Annotation DMLs after filtering to keep only standard chromosomes (",nrow(annotated_dml),"loci)\n")
# Annotation DMLs after filtering to keep only standard chromosomes ( 202 loci)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Step 6: GO Enrichment Analysis and KEGG pathway
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract unique Entrez Gene IDs
dml_genes = unique(annotated_dml$geneId)
dml_genes = dml_genes[!is.na(dml_genes)] # Remove NAs

len(dml_genes)

# -----------------------------------------------------------
# Step 6.1: GO enrichment analysis 
# ----------------------------------------------------------
  # GO enrichment analysis (Biological Process)
go_bp_results = enrichGO(dml_genes, 
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)

# Extract Enriched terms (if any exist)
# Convert enrichment results to dataframe
go_bp_df = as.data.frame(go_bp_results)

# Check if there are any enriched terms
if (nrow(go_bp_df) > 0) {
  
  # Print the top 10 enriched GO terms
  print(head(go_bp_df, 10))
} else {
  cat ("No enriched GO terms (BP) found.")
}

# No enriched GO terms (BP) found.

# GO enrichment analysis (Molecular Functions) 
go_mf_results = enrichGO(dml_genes,
                         OrgDb = org.Hs.eg.db,
                         keyType =  "ENTREZID",
                         ont = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)

# Extract enriched terms (if any exist)
go_mf_df = as.data.frame(go_mf_results)

# Check if there are any enriched terms
if (nrow(go_mf_df) > 0) {
  
  # Print the top 5 enriched GO terms
  print (head(go_mf_df, 10))
} else {
  cat ("No enriched GO terms (MF) found.")
}
# No enriched GO terms (MF) found.

# ---------------------------------------------------------------
# Step 6.2: KEGG Pathway
# -------------------------------------------------------------------
# KEGG PATHWAY ENRICHMENT
kegg_results <- enrichKEGG(gene = dml_genes, 
                           organism = "hsa",  # Homo sapiens, 
                           pvalueCutoff = 0.05)

# Extract KEGG enrichment terms
kegg_df = as.data.frame(kegg_results)

# Check if there are any enriched terms
if (nrow(kegg_df) > 0) {
  # Print the top 5 enriched KEGG terms
  print(head(kegg_df, 5))
} else {
  cat ("No enriched KEGG terms")
}

#                   category            subcategory    ID                   Description       GeneRatio BgRatio RichFactor FoldEnrichment   zScore       pvalue   p.adjust     qvalue               geneID Count
# hsa05217     Human Diseases Cancer: specific types hsa05217           Basal cell carcinoma      4/53  63/8541 0.06349206       10.23181 5.811315 0.0006004906 0.04703305 0.04462735 1855/7481/83439/6608     4
# hsa04966 Organismal Systems       Excretory system hsa04966 Collecting duct acid secretion      3/53  28/8541 0.10714286       17.26617 6.812198 0.0006624373 0.04703305 0.04462735       10312/496/6521     3

# =====================================================================================================
# Let's delve into the genes associated with the two enriched pathways identified in your KEGG pathway enrichment analysis: Basal cell carcinoma (hsa05217) and Collecting duct acid secretion (hsa04966).

# 1. Basal Cell Carcinoma (hsa05217):
# Basal cell carcinoma (BCC) is a common form of skin cancer primarily linked to aberrant activation of the Hedgehog signaling pathway.
# Mutations in key genes such as PTCH1 (Gene ID: 5727), SMO  (Gene ID: 6608), and SHH (Gene ID: 6469) lead to continuous activation of target genes, promoting uncontrolled cell proliferation. Additionally, mutations in the tumor suppressor gene TP53 are frequently observed in BCC cases. 

# In our analysis, the following genes were implicated in this pathway:

# PTCH1 (Gene ID: 5727): Acts as a receptor in the Hedgehog signaling pathway, regulating cell growth.

# SMO (Gene ID: 6608): A G protein-coupled receptor that transduces Hedgehog signals, influencing cell proliferation.

# SHH (Gene ID: 6469): Encodes the Sonic Hedgehog protein, a key signaling molecule in developmental processes.

# TP53 (Gene ID: 7157): A tumor suppressor gene involved in DNA repair and apoptosis; its mutation is common in various cancers, including BCC.


# 2. Collecting Duct Acid Secretion (hsa04966):
# The collecting duct of the kidney plays a crucial role in maintaining acid-base balance by secreting hydrogen ions. 
# This process involves various transporters and enzymes that regulate urine acidity. 


# The genes from our analysis associated with this pathway include:

# ATP6V1B1 (Gene ID: 10312): Encodes a subunit of the vacuolar ATPase involved in proton transport.

# SLC4A1 (Gene ID: 6521): Encodes an anion exchanger critical for bicarbonate transport in renal cells.

# CA2 (Gene ID: 760): Encodes carbonic anhydrase II, an enzyme that catalyzes the reversible hydration of carbon dioxide, essential in acid-base balance.

# These genes are integral to the physiological processes of their respective pathways. 
# Alterations or dysregulation in these genes can contribute to disease development, such as cancer progression or renal dysfunction.
# ===============================================================================================================================================================

# Save KEGG results
write.csv(kegg_df, "KEGG_enrichment_DMLs.csv", row.names = FALSE)

# Visualize KEGG enrichment results
kegg_dotplot = dotplot(kegg_results, showCategory = nrow(kegg_results), title = "KEGG Pathway Enrichment of DMLs")

# Save the figure
ggsave("KEGG_Pathway_Enrichment_DMLs.png", plot = kegg_dotplot, height = 6, dpi = 300)

# --------------------------------------------------------------------------------------------------------------------------
# Step 7: Visualization
# ----------------------------------------------------------------------------------------------------------------------------

# 7.1 **Volcano plots of DMLs**
ggplot(dml_test, 
       aes(x = stat, 
           y = -log10(pvals),
           color = fdrs < 0.05)) +
  geom_point (alpha = 0.7) +
  scale_color_manual(values = c("black", "red")) +
  labs (title = "Volcano Plots of DMLs",
        x = "Methylation Difference (Stat)",
        y = "-log10(p-value)") +
  theme_minimal()
ggsave("Volcano_plot_DMLs.png")

# 7.2 **Heatmap of Significant DMLs**
# Extract per-sample methylation values for significant DMLs
meth_values <- getMeth(BSobj_filtered, type = "raw")

# Verify the dimensions of meth_values
cat("Methylation matrix dimensions:", dim(meth_values), "\n")
# Methylation matrix dimensions: 1874421 6 

# Ensure row names in `significant_dml` match chromosome-position format
significant_dml$site <- paste0(significant_dml$chr, ":", significant_dml$pos)
bsobj_sites <- paste0(seqnames(rowRanges(BSobj_filtered)), ":", start(rowRanges(BSobj_filtered)))

# Set row names in meth_values to match CpG site identifiers
rownames(meth_values) <- bsobj_sites


if (nrow(significant_dml) > 1) {
  common_sites = intersect(rownames(meth_values), significant_dml$site)
  filtered_meth_values = meth_values[common_sites, ]
  
  # Generate heatmap
  pheatmap(filtered_meth_values, cluster_rows = TRUE, cluster_cols = TRUE,
           color = colorRampPalette(c("blue", "white", "red")) (20),
           main = "Heatmap of Significant DMLs")
  ggsave("heatmap_DMLs_across_samples.png", width = 8, height = 6, dpi = 300)
}
  
# We can clearly see that pre-chemo samples (SRR19127197, SRR19127202 and SRR191127198) are very similar.
# Post-chemo samples (SRR19127182 and SRR19127178) are very similar while SRR19127177 is less similar. 

# -------------------------------------------------------------------------------------------------------------------
# B. DMR Analysis with Pathway Enrichment
# --------------------------------------------------------------------------------------------------------------------  

# =======================================================================================================
# DMR ANALYSIS PIPELINE - DSS & Functional Annotation
# =======================================================================================================
# This script performs Differentially Methylated Region (DMR) analysis using the DSS package in R.
# The workflow includes:
# 1. Data Preprocessing
# 2. Differential Methylation Region (DMR) Analysis
# 3. Functional Annotation of DMRs (Gene Mapping)
# 4. Pathway & Gene Ontology (GO) Enrichment
# 5. Visualization of DMRs
# =======================================================================================================

# -------------------------------------------------------------------------------------------------------
# INSTALL REQUIRED PACKAGES (If not already installed)
# -------------------------------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DSS", "bsseq", "ggplot2", "reshape2", "pheatmap", "GenomicRanges",
                       "biomaRt", "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", 
                       "org.Hs.eg.db", "clusterProfiler", "Gviz"), ask = FALSE)

# -------------------------------------------------------------------------------------------------------
# LOAD NECESSARY LIBRARIES
# -------------------------------------------------------------------------------------------------------
library(DSS)           # DMR analysis
library(bsseq)         # Methylation data processing
library(ggplot2)       # Visualization
library(reshape2)      # Data reshaping
library(pheatmap)      # Heatmaps
library(GenomicRanges) # Genomic operations
library(biomaRt)       # Gene annotation
library(ChIPseeker)    # Peak annotation
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Human genome annotation
library(org.Hs.eg.db)  # Gene symbol mapping
library(clusterProfiler) # Functional enrichment analysis
library(Gviz)          # Genome visualization

# -------------------------------------------------------------------------------------------------------
# SET WORKING DIRECTORY (Modify this path accordingly)
# -------------------------------------------------------------------------------------------------------
setwd("D:/My projects_Oct 2023/Bisulfite_seq/Bisulfite Sequencing")

# =======================================================================================================
# 1. LOAD & PREPROCESS DATA
# =======================================================================================================

# Load DSS input file (methylation counts for CpG sites)
dss_data <- read.csv("D:/My projects_Oct 2023/Bisulfite_seq/Bisulfite Sequencing/dss_input.csv")  # Generated from Python preprocessing

# Ensure numeric values for methylation counts
dss_data$Total_C <- as.numeric(as.character(dss_data$Total_C))
dss_data$Methylated_C <- as.numeric(as.character(dss_data$Methylated_C))

# Convert data to BSseq format
sample_list <- split(dss_data, dss_data$Sample_ID)

bsseq_input <- lapply(sample_list, function(df) {
  data.frame(chr = df$Chromosome, pos = df$Position, N = df$Total_C, X = df$Methylated_C)
})

sample_names <- unique(dss_data$Sample_ID)

# Create BSseq object
BSobj <- makeBSseqData(bsseq_input, sampleNames = sample_names)

# =======================================================================================================
# 2. ASSIGN SAMPLE METADATA & FILTER LOW-COVERAGE SITES
# =======================================================================================================

# Create metadata (design matrix)
design_matrix <- data.frame(
  Sample_ID = sampleNames(BSobj),
  Condition = factor(c("Pre_chemo", "Post_chemo", "Pre_chemo", "Pre_chemo", "Post_chemo", "Post_chemo"))
)
rownames(design_matrix) <- sampleNames(BSobj)
pData(BSobj) <- design_matrix

# Filter low-coverage sites (CpG sites with at least 3 reads)
BSobj_filtered <- BSobj[rowSums(getCoverage(BSobj)) > 3, ]  

# Check for missing values & remove NAs
#BSobj_filtered <- BSobj_filtered[complete.cases(getCoverage(BSobj_filtered)), ]

nrow(BSobj_filtered)

# =======================================================================================================
# 3. PERFORM DIFFERENTIAL METHYLATION REGION (DMR) ANALYSIS
# =======================================================================================================

# Identify differentially methylated regions (DMRs)
dmr_results <- callDMR(BSobj_filtered, p.threshold = 0.05)

# Save DMR results
write.csv(dmr_results, "significant_DMRs.csv", row.names = FALSE)
cat("Significant DMRs saved (", nrow(dmr_results), " regions)\n")# 14

#dmr_results = read.csv("significant_DMRs.csv")
# =======================================================================================================
# 4. FUNCTIONAL ANNOTATION OF DMRs
# =======================================================================================================

# Load genome annotation
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Convert DMRs to GRanges format
gr_dmr <- GRanges(seqnames = dmr_results$chr,
                  ranges = IRanges(start = dmr_results$start, 
                                   end = dmr_results$end))

# Annotate DMRs (Assign each region to its nearest gene, promoter, exon, or intron)
peak_annotation <- annotatePeak(gr_dmr, 
                                TxDb = txdb, 
                                tssRegion = c(-2000, 2000),
                                verbose = FALSE)

# Convert to data frame
annotated_dmr <- as.data.frame(peak_annotation)

# Save annotated results
write.csv(annotated_dmr, "annotated_DMRs.csv", row.names = FALSE)
cat("Annotated DMRs saved as 'annotated_DMRs.csv'\n")

# =======================================================================================================
# 5. PATHWAY & GENE ONTOLOGY (GO) ENRICHMENT ANALYSIS
# =======================================================================================================

# Extract unique Entrez Gene IDs
dmr_genes <- unique(annotated_dmr$geneId)
dmr_genes <- dmr_genes[!is.na(dmr_genes)]  # Remove NAs

# Perform GO Enrichment Analysis (Biological Process)
go_bp_results <- enrichGO(gene = dmr_genes,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

# Save GO results
write.csv(as.data.frame(go_bp_results), "GO_Enrichment_BP_DMRs.csv", row.names = FALSE)
cat("GO Enrichment results saved as 'GO_Enrichment_DMRs.csv'\n")

# Perform GO Enrichment Analysis (molecular functions)
go_mf_results <- enrichGO(gene = dmr_genes,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "MF",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
# Save GO results
write.csv(as.data.frame(go_mf_results), "GO_Enrichment_MF_DMRs.csv", row.names = FALSE)
cat("GO EMnrichment results saved as 'GO_Enrichment_DMRs.csv'\n")

# Perform KEGG Pathway Analysis
kegg_results <- enrichKEGG(gene = dmr_genes,
                           organism = "hsa",
                           pvalueCutoff = 0.05)
# 0 enriched terms found
# Save KEGG results
write.csv(as.data.frame(kegg_results), "KEGG_Enrichment_DMRs.csv", row.names = FALSE)
cat("KEGG Enrichment results saved as 'KEGG_Enrichment_DMRs.csv'\n")

# =======================================================================================================
# 6. VISUALIZATION OF DMRs
# =======================================================================================================

# Define genome assembly
genome <- "hg38"

# Select a chromosome for visualization
selected_chr <- "chr7"

# Filter DMRs for the selected chromosome
gr_dmr_filtered <- gr_dmr[as.character(seqnames(gr_dmr)) == selected_chr]

# Ensure DMRs exist for the selected chromosome
if (length(gr_dmr_filtered) == 0) {
  stop(paste("⚠ No DMRs found on", selected_chr, "- Try another chromosome."))
}

# Assign dummy score for visualization
mcols(gr_dmr_filtered)$score <- rep(1, length(gr_dmr_filtered))

# Create DMR Track
dmr_track <- DataTrack(range = gr_dmr_filtered, genome = genome, chromosome = selected_chr,
                       name = "DMRs", type = "histogram", col.histogram = "blue", 
                       fill.histogram = "lightblue", lwd = 2)

# Create Genome Axis Track
axis_track <- GenomeAxisTrack()

# Define Plot Range
plot_start <- min(start(gr_dmr_filtered)) - 50  
plot_end <- max(end(gr_dmr_filtered)) + 50  

# Generate Plot
plotTracks(list(axis_track, dmr_track), 
           from = plot_start, to = plot_end, 
           col.title = "black", background.title = "lightgray",
           main = paste("Differentially Methylated Regions on", selected_chr))

# =======================================================================================================
# SUMMARY OF KEY OUTPUTS
# =======================================================================================================
cat("\n **Summary of DMR Analysis:**\n")
cat(" Number of significant DMRs: ", nrow(dmr_results), "\n")
cat(" Functional annotation saved: 'annotated_DMRs.csv'\n")
cat(" GO & KEGG enrichment results saved\n")
cat(" Visualization completed\n")

# ------------------------END OF GO & KEGG ENRICHMENT ANALYSIS ON DML and DMR -------------------------------------------------------------------

  
  
  
  
  
  
  
  
  
  
  