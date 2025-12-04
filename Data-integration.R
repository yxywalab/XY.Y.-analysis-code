
#' Integrate multiple single-cell RNA-seq datasets and analyze cell type composition
#' 
#' This function performs integration of multiple single-cell datasets using Seurat's
#' integration pipeline, followed by clustering, visualization, and compositional analysis
#' of cell types across different data sources.
#' 
#' @details
#' The pipeline consists of four main steps:
#' 1. Parallel processing setup for efficient computation
#' 2. Dataset integration using reciprocal PCA (RPCA) and anchor-based integration
#' 3. Clustering, dimensionality reduction, and visualization
#' 4. Cell type composition analysis across different data sources
#' 
#' @note
#' This script requires substantial computational resources (up to 500GB RAM).
#' Adjust parameters according to available system resources.

# Enable parallel processing to speed up computations
# Sets up multicore parallelization with 10 worker processes
library(future)
plan("multicore", workers = 10) 

# Increase maximum allowed size for global variables in parallel processing
# Allows transfer of large objects (up to 500GB) between parallel processes
options(future.globals.maxSize = 500 * 1024^3) 

# Create list of Seurat objects for integration
# Combines multiple datasets (brain and other tissues) into a single list for batch correction
seurat_list <- list(scRNA_brain, scRNA1, scRNA3, scRNA4) 

#' Identify integration anchors between datasets
#' Uses reciprocal PCA (RPCA) to find mutual nearest neighbors between datasets
#' RPCA is more robust to batch effects than CCA for large dataset integration
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,      # List of Seurat objects to integrate
  anchor.features = 2000,         # Number of variable features to use for integration
  dims = 1:30,                    # Dimensions from PCA to use for anchor finding
  k.anchor = 20,                  # Number of neighbors to use in anchor filtering
  reduction = "rpca",             # Use reciprocal PCA for anchor finding
  verbose = FALSE                 # Suppress progress messages
)

# Pause execution to ensure proper memory cleanup
# Allows garbage collection and prevents memory conflicts
Sys.sleep(2)

# Remove intermediate objects to free up memory
# These objects are no longer needed after anchor identification
rm(scRNA_all, scRNA_brain_epi, scRNA_tlung)

#' Integrate datasets using identified anchors
#' Creates a combined assay with batch-corrected expression values
#' This step removes technical batch effects while preserving biological variation
integrated <- IntegrateData(
  anchorset = anchors,            # Anchor set from FindIntegrationAnchors
  dims = 1:30,                    # Dimensions to use for integration
  verbose = FALSE                 # Suppress progress messages
)

# Additional pause for memory management
Sys.sleep(2)

# Set the integrated assay as default for downstream analysis
# All subsequent operations will use batch-corrected expression values
DefaultAssay(integrated) <- "integrated"

#' Scale and center the integrated data
#' Normalizes expression values for each gene across all cells
#' Required for dimensionality reduction techniques
integrated <- ScaleData(
  integrated,
  features = VariableFeatures(integrated))  # Use variable features identified during integration

#' Perform Principal Component Analysis (PCA)
#' Reduces dimensionality of the integrated data for downstream analysis
integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated))

#' Run Uniform Manifold Approximation and Projection (UMAP)
#' Non-linear dimensionality reduction for visualization
# Uses fewer neighbors (n.neighbors=5) for local structure preservation
integrated <- RunUMAP(integrated, n.neighbors = 5, dims = 1:20)

# Visualize PCA results
PCAPlot(integrated) 

#' Construct shared nearest neighbor (SNN) graph
#' Builds k-nearest neighbor graph using first 10 principal components
# Lower dimensionality (dims=1:10) focuses on major sources of variation
integrated <- FindNeighbors(integrated, dims = 1:10) 

#' Perform graph-based clustering
#' Identifies cell communities in the SNN graph using Louvain algorithm
# Low resolution (0.1) produces broad cell type clusters
integrated <- FindClusters(integrated, dims = 1:10, resolution = c(0.1))

#' Generate UMAP visualization with optimized parameters
# Fewer dimensions (1:5) for clearer separation of major clusters
integrated <- RunUMAP(integrated, dims = 1:5)

# Visualize clustering results with custom colors
DimPlot(integrated, reduction = "umap") +
  scale_color_manual(values = colors) +
  plot_annotation(title = "Integrated Dataset Clustering")

#' Load pre-processed Seurat object for detailed cell type analysis
# Contains annotated cell types from previous analysis
scRNA <- readRDS("ALL_seurat.rds")  

# Define color palette for cell type visualization
# Using Set2 palette from RColorBrewer for colorblind-friendly visualization
my_colors <- RColorBrewer::brewer.pal(n = 8, name = "Set2")
# Color codes: "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"

# Visualize cell types in UMAP space
DimPlot(scRNA, 
        reduction = 'umap',                 # Use UMAP coordinates
        group.by = 'Celltype',             # Color by cell type annotation
        pt.size = 1,                       # Point size for cells
        label = T,                         # Add cluster labels
        label.size = 5,                    # Label text size
        alpha = c(0.6)) +                  # Transparency for overlapping points
  theme(legend.position = "none") +        # Hide legend (labels provide identification)
  scale_color_manual(values = my_colors)   # Apply custom color palette

# Check unique cell types present in the dataset
unique(scRNA$celltype)

# Convert celltype to character for string manipulation
scRNA$celltype <- as.character(scRNA$celltype)

#' Clean and standardize cell type annotations
# Reclassify "Smooth_muscle_cells" and NA values as "Fibroblast"
scRNA$celltype[is.na(scRNA$celltype) | scRNA$celltype == "Smooth_muscle_cells"] <- "Fibroblast"

# Convert back to factor for efficient storage and plotting
scRNA$celltype <- factor(scRNA$celltype)

# Check unique sample identifiers
unique(scRNA$orig.ident)

# Load additional required libraries
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggpubr)

#' Define data source mapping
# Maps sample identifiers (orig.ident) to broader data source categories
# Each source represents a distinct study or data generation batch
sources <- list(
  source1 = "SeuratProject",  # Single sample
  source2 = c("GSM7475327", "GSM7475328"),  # Two samples from same study
  source3 = c("GSM5645894", "GSM5645895", "GSM5645896"),  # Three samples
  source4 = c("GSM6112137_LC1", "GSM6112138_LC2", "GSM6112139_LC3", 
              "GSM6112140_LC4", "GSM6112141_LC5", "GSM6112142_LC6",
              "GSM6112143_LC7", "GSM6112144_LC8", "GSM6112145_LC9", 
              "GSM6112146_LC10"),  # 10 lung cancer samples
  source5 = c(  # 29 NSCLC samples from different patients
    "GSM6957590_NSCLC_PA001_sn_raw_feature_bc_matrix",
    "GSM6957591_NSCLC_PA004_sn_raw_feature_bc_matrix",
    "GSM6957592_NSCLC_PA005_sn_raw_feature_bc_matrix",
    # ... (additional samples)
    "GSM6957628_NSCLC_KRAS_17_sn_raw_feature_bc_matrix"
  ),
  source6 = c(  # 41 samples from OMIX dataset
    "OMIX007088-01", "OMIX007088-03", "OMIX007088-07", "OMIX007088-11",
    # ... (additional samples)
    "OMIX007088-85"
  )
)

# Initialize source column with missing values
scRNA$source <- NA_character_

# Assign source categories based on sample identifiers
# Loops through each source definition and matches orig.ident values
for (i in seq_along(sources)) {
  sample_names <- sources[[i]]
  scRNA$source[scRNA$orig.ident %in% sample_names] <- names(sources)[i]
}

# Check for samples that weren't assigned to any source
# Useful for debugging mapping completeness
unassigned <- scRNA@meta.data %>% 
  filter(is.na(source)) %>% 
  distinct(orig.ident)

#' Calculate cell type proportions by data source
# Counts cells for each source-celltype combination and calculates percentages
cell_prop <- scRNA@meta.data %>%
  group_by(source, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(source) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Load additional libraries for advanced plotting
library(ggplot2)
library(scales)
library(dplyr)

#' Create stacked bar plot of cell type composition
# Visualizes proportional distribution of cell types across different data sources
p_stacked <- ggplot(cell_prop, aes(x = source, y = prop, fill = celltype)) +
  
  # Create stacked bars with controlled width and reverse stacking order
  geom_col(width = 0.65, position = position_stack(reverse = TRUE)) +
  
  # Add percentage labels for proportions >5%
  # Position labels in the middle of each segment
  geom_text(aes(label = ifelse(prop > 0.05, percent(prop, accuracy = 1), "")), 
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 3.5, color = "black", show.legend = FALSE) +

  # Configure y-axis with percentage formatting
  scale_y_continuous(expand = c(0, 0),  # Remove padding at axis limits
                     labels = percent_format(accuracy = 1),  # Format as percentages
                     breaks = seq(0, 1, 0.2),  # Major gridlines every 20%
                     limits = c(0, 1.05)) +  # Slight extension for label spacing
  
  # Add padding to x-axis for better visual balance
  scale_x_discrete(expand = expansion(add = 0.5)) +
  
  # Apply custom color scheme for cell types
  scale_fill_manual(values = custom_colors) +
  
  # Customize plot theme for publication-quality output
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(color = "black"),  # All text in black
    
    # Axis line and tick customization
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    
    # Axis title adjustments
    axis.title.x = element_blank(),  # Remove x-axis title (source names are self-explanatory)
    axis.title.y = element_text(size = 15, margin = margin(r = 10), vjust = 2),
    
    # X-axis text with 45-degree angle for readability
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,
      vjust = 1,
      size = 15,
      color = "black"
    ),
    
    # Y-axis text styling
    axis.text.y = element_text(size = 15, color = "black"),
    
    # Legend positioning and styling
    legend.position = "right",
    legend.title = element_blank(),  # Remove legend title
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.text = element_text(size = 15),
    legend.background = element_rect(fill = "white", color = NA),
    
    # Plot margins and background
    plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"),  # Top, right, bottom, left
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),  # Remove gridlines for cleaner look
    plot.title = element_text(hjust = 0.5, size = 14)  # Centered title
  ) +
  
  # Plot titles and labels
  labs(
    title = "Cell Type Composition by Data Source",
    y = "Percentage of Cells"
  )

# Display the final plot
print(p_stacked)
