#' Calculate cell type diversity and distribution similarity metrics
#'
#' This function analyzes single-cell RNA-seq metadata to calculate:
#' 1. Shannon diversity index for each data source
#' 2. Jensen-Shannon divergence between each source and average distribution
#' 3. Multidimensional scaling (MDS) of cell type distributions
#'
#' @param scRNA A Seurat object containing single-cell RNA-seq data with
#' metadata columns 'source' (data source identifier) and 'celltype' (cell type annotation)
#'
#' @return A list containing three ggplot objects:
#' \itemize{
#' \item{p_diversity}{Bar plot of Shannon diversity indices by source}
#' \item{p_similarity}{Bar plot of similarity scores to average distribution}
#' \item{p_mds}{MDS plot showing distribution of sources in 2D space}
#' }
#'
#' @details
#' The function performs three complementary analyses:
#' 1. Diversity Analysis: Calculates Shannon diversity index to quantify
#' the richness and evenness of cell types within each data source
#' 2. Similarity Analysis: Computes Jensen-Shannon divergence between each
#' source's cell type distribution and the average distribution across all sources
#' 3. Dimensionality Reduction: Performs MDS to visualize sources in 2D space
#' based on Euclidean distances of their cell type distributions
#'
#' @examples
#' \dontrun{
#' results <- analyze_celltype_diversity(scRNA_object)
#' results$p_diversity | results$p_similarity | results$p_mds
#' }
#'
#' @export
library(vegan) 
library(vegan)    
library(philentropy) 
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggrepel)
cell_prop <- scRNA@meta.data %>% dplyr::count(source, celltype) %>% dplyr::group_by(source) %>%dplyr::mutate(prop = n / sum(n)) %>%dplyr::ungroup()
diversity_data <- cell_prop %>%
  tidyr::pivot_wider(id_cols = source,names_from = celltype, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("source") %>% vegan::diversity(index = "shannon") %>% as.data.frame() %>%
  tibble::rownames_to_column("source") %>% dplyr::rename(shannon = 2) %>% dplyr::arrange(-shannon)
p_diversity <- ggplot(diversity_data, aes(x = reorder(source, shannon), y = shannon)) +
  geom_bar(stat = "identity", fill = "#2a9d8f", width = 0.6) +
  geom_text(aes(label = round(shannon, 2)), vjust = -0.5, size = 5, fontface = "bold") +
  geom_hline(yintercept = max(diversity_data$shannon), linetype = "dashed", color = "#e76f51", size = 1) +
  annotate("text", x = 3, y = max(diversity_data$shannon) + 0.1, 
           label = paste("Highest Diversity:"), 
           color = "#e76f51", size = 5, fontface = "bold") +
  labs(title = "Cell Type Diversity by Source",  
       x = "Data Source", y = "Shannon Diversity Index") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        plot.title = element_text(face = "bold"))
avg_distribution <- cell_prop %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(avg_prop = mean(prop)) %>%
  dplyr::ungroup()
jsd_matrix <- cell_prop %>%
  tidyr::pivot_wider(id_cols = source, 
                     names_from = celltype, 
                     values_from = prop, 
                     values_fill = 0) %>%
  tibble::column_to_rownames("source") %>%
  philentropy::distance(method = "jensen-shannon") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("source") %>%
  dplyr::rename(jsd = 2) %>%
  dplyr::mutate(
    similarity = 1 - jsd,
    rank = dense_rank(jsd)
  )
p_similarity <- ggplot(jsd_matrix, aes(x = reorder(source, -similarity), y = similarity)) +
  geom_bar(stat = "identity", fill = "#7b2cbf", width = 0.6) +
  geom_text(aes(label = paste0(round(similarity, 3), "\nRank ", rank, "")), 
            vjust = 1.2, size = 4, fontface = "bold") +
  geom_hline(yintercept = max(jsd_matrix$similarity), linetype = "dashed", color = "#e76f51", size = 1) +
  labs(title = "Similarity to Average Cell Type Distribution",
       x = "Data Source", y = "Similarity Score") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        plot.title = element_text(face = "bold"))
mds_data <- cell_prop %>%
  tidyr::pivot_wider(id_cols = source, 
                     names_from = celltype, 
                     values_from = prop, 
                     values_fill = 0) %>%
  tibble::column_to_rownames("source")
dist_matrix <- dist(mds_data, method = "euclidean")
mds_result <- cmdscale(dist_matrix, k = 2, eig = TRUE) %>%
  .$points %>%
  as.data.frame() %>%
  tibble::rownames_to_column("source") %>%
  dplyr::rename(MDS1 = V1, MDS2 = V2)
centroid <- mds_result %>%
  dplyr::summarise(MDS1 = mean(MDS1), MDS2 = mean(MDS2))
mds_result <- mds_result %>%
  dplyr::mutate(
    dist_to_center = sqrt((MDS1 - centroid$MDS1)^2 + (MDS2 - centroid$MDS2)^2),
    centrality_rank = dense_rank(dist_to_center)
  )
p_mds <- ggplot(mds_result, aes(x = MDS1, y = MDS2, color = dist_to_center)) +
  geom_point(size = 8) +
  geom_point(data = centroid, aes(x = MDS1, y = MDS2), 
             shape = 8, size = 10, color = "#e76f51", stroke = 1.5) +
  geom_text_repel(aes(label = paste0(source, "\n(Rank ", centrality_rank, ")")), 
                  size = 5, fontface = "bold", box.padding = 0.8) +
  geom_segment(aes(xend = centroid$MDS1, yend = centroid$MDS2), 
               linetype = "dashed", alpha = 0.6) +
  annotate("text", x = min(mds_result$MDS1)+0.3, y = max(mds_result$MDS2), 
           label = "Closer to center = More representative", 
           hjust = 0, size = 5, color = "#e76f51", fontface = "bold") +
  scale_color_gradient(low = "#2a9d8f", high = "#e76f51", name = "Distance to Center") +
  labs(title = "Multidimensional Scaling of Cell Type Distributions") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

