#' Calculate EMT signatures using AUCell and AddModuleScore methods
#' 
#' This function computes Epithelial-Mesenchymal Transition (EMT) signatures
#' for single-cell RNA-seq data using two complementary methods: AUCell and
#' Seurat's AddModuleScore. It integrates gene sets from MSigDB and provides
#' visualization of EMT activity across cells.
#' 
#' @details
#' The function performs the following steps:
#' 1. Retrieves EMT gene signatures from MSigDB (Hallmark and GO categories)
#' 2. Computes EMT activity scores using AUCell method
#' 3. Computes EMT module scores using Seurat's AddModuleScore
#' 4. Visualizes the results using FeaturePlots and Violin plots
#' 
#' @note
#' The AUCell method ranks cells based on expression of gene sets and calculates
#' the Area Under the Curve (AUC) for each cell. AddModuleScore calculates the
#' average expression of a gene set subtracted by the average expression of
#' control gene sets.

# Load required libraries for single-cell analysis and gene set enrichment
library(Seurat)      # For single-cell RNA-seq analysis
library(msigdbr)     # For accessing MSigDB gene sets
library(AUCell)      # For gene set enrichment analysis at single-cell level
library(dplyr)       # For data manipulation
library(patchwork)   # For combining multiple plots

#' Retrieve EMT gene signatures from MSigDB
#' MSigDB (Molecular Signatures Database) contains curated gene sets

# Get Hallmark gene sets for human
msig_h <- msigdbr(species = "Homo sapiens", category = "H")

# Extract EMT signature from Hallmark collection
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION is a curated set of 200 genes
emt_genes_hallmark <- msig_h %>%
  filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>%
  pull(gene_symbol)

# Print the list of Hallmark EMT genes for verification
print(emt_genes_hallmark)

# Get GO (Gene Ontology) gene sets for human
msig_c5 <- msigdbr(species = "Homo sapiens", category = "C5")

# Extract EMT signature from GO Biological Processes
# GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION contains genes involved in EMT process
emt_genes_go <- msig_c5 %>%
  filter(gs_name == "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION") %>%
  pull(gene_symbol)

#' Combine and prepare EMT gene list
# Merge Hallmark and GO EMT gene sets, remove duplicates
emt_genes <- unique(c(emt_genes_hallmark, emt_genes_go)) %>% 
  toupper()  # Convert to uppercase for consistency with gene names in scRNA object

#' Prepare expression matrix for AUCell analysis
# Extract raw count matrix from Seurat object
# AUCell requires a matrix with genes as rows and cells as columns
expr_matrix <- as.matrix(GetAssayData(scRNA, slot = "counts"))

#' Perform AUCell analysis
# AUCell calculates whether critical subset of genes in the gene set are 
# enriched in the top expressed genes for each cell

# Step 1: Build rankings of genes for each cell
# Cells are ranked by expression level, highest expressed genes get top ranks
cells_rankings <- AUCell_buildRankings(expr_matrix, nCores = 4)  

# Step 2: Calculate AUC scores for EMT gene set
# aucMaxRank: Only consider top 5% of expressed genes for AUC calculation
# This focuses on highly expressed genes which are most biologically relevant
cells_AUC <- AUCell_calcAUC(geneSets = list(EMT = emt_genes), 
                            rankings = cells_rankings, 
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

# Add AUCell scores to Seurat object metadata
# Extract AUC scores and add as a new metadata column
scRNA5$EMT_AUCell <- as.numeric(getAUC(cells_AUC)["EMT", ])

#' Prepare gene list for AddModuleScore
# Filter EMT genes to only include those present in the scRNA dataset
emt_genes_filtered <- intersect(emt_genes, rownames(scRNA))

#' Calculate AddModuleScore for EMT signature
# AddModuleScore calculates the average expression of a set of genes, 
# subtracted by the average expression of control feature sets
# ctrl = 100: Use 100 control gene sets of similar expression levels
scRNA <- AddModuleScore(scRNA,
                        features = list(EMT = emt_genes_filtered),
                        name = "EMT_Module",
                        ctrl = 100)

# Rename the module score column for clarity
# The function creates a column named EMT_Module1 by default
colnames(scRNA@meta.data) <- gsub("EMT_Module1", "EMT_AddModuleScore", 
                                  colnames(scRNA@meta.data))

#' Visualize EMT scores using FeaturePlot
# Create side-by-side UMAP plots showing both scoring methods
# order = TRUE: Plot cells with highest scores on top for better visualization
FeaturePlot(scRNA, features = c("EMT_AUCell", "EMT_AddModuleScore"), 
            order = TRUE, combine = FALSE, reduction = "umap") %>%
  wrap_plots(ncol = 2)

# Alternative visualization with custom color palette
# Blue (low) -> Yellow (medium) -> Red (high) gradient
plot_list <- FeaturePlot(scRNA,
                         features = c("EMT_AUCell", "EMT_AddModuleScore"), 
                         order = TRUE,
                         combine = FALSE,
                         reduction = "umap",
                         cols = c("blue", "yellow", "red")) 

# Combine the two plots side by side
wrap_plots(plot_list, ncol = 2)

#' Compare EMT scores across clusters using Violin plot
# Visualize distribution of AUCell scores across Louvain clusters
# pt.size = 0: Remove individual points for cleaner visualization
VlnPlot(scRNA, features = "EMT_AUCell", group.by = "RNA_snn_res.0.2", pt.size = 0)
