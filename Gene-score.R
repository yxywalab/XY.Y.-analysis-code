library(Seurat)
library(msigdbr)
library(AUCell)
library(dplyr)
library(patchwork)

msig_h <- msigdbr(species = "Homo sapiens", category = "H")
emt_genes_hallmark <- msig_h %>%
  filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>%
  pull(gene_symbol)
print(emt_genes_hallmark )

msig_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
emt_genes_go <- msig_c5 %>%
  filter(gs_name == "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION") %>%
  pull(gene_symbol)

emt_genes <- unique(c(emt_genes_hallmark, emt_genes_go)) %>% 
  toupper() 

expr_matrix <- as.matrix(GetAssayData(scRNA, slot = "counts")) 

cells_rankings <- AUCell_buildRankings(expr_matrix, nCores = 4)  
cells_AUC <- AUCell_calcAUC(geneSets = list(EMT = emt_genes), 
                            rankings = cells_rankings, 
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))  

scRNA5$EMT_AUCell <- as.numeric(getAUC(cells_AUC)["EMT", ])

emt_genes_filtered <- intersect(emt_genes, rownames(scRNA))

scRNA <- AddModuleScore(scRNA,
                             features = list(EMT = emt_genes_filtered),
                             name = "EMT_Module",
                             ctrl = 100) 

colnames(scRNA@meta.data) <- gsub("EMT_Module1", "EMT_AddModuleScore", 
                                       colnames(scRNA_15@meta.data))

FeaturePlot(scRNA, features = c("EMT_AUCell", "EMT_AddModuleScore"), 
            order = TRUE, combine = FALSE,reduction = "umap") %>%
  wrap_plots(ncol = 2)
plot_list <- FeaturePlot(scRNA,
                         features = c("EMT_AUCell", "EMT_AddModuleScore"), 
                         order = TRUE,
                         combine = FALSE,
                         reduction = "umap",
                         cols = c("blue", "yellow", "red")) 

wrap_plots(plot_list, ncol = 2)

VlnPlot(scRNA, features = "EMT_AUCell", group.by = "RNA_snn_res.0.2", pt.size = 0)
