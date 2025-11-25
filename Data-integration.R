
library(future)
plan("multicore", workers = 10) 

options(future.globals.maxSize = 500 * 1024^3) 

seurat_list <- list(scRNA_brain, scRNA1, scRNA3,scRNA4) 

anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = 2000,
  dims = 1:30,
  k.anchor = 20, 
  reduction = "rpca", 
  verbose = FALSE
)
Sys.sleep(2)
rm(scRNA_all,scRNA_brain_epi,scRNA_tlung)

integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:30,
  verbose = FALSE
)
Sys.sleep(2)
DefaultAssay(integrated) <- "integrated"

integrated <- ScaleData(
  integrated,
  features = VariableFeatures(integrated))
integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated))
integrated <- RunUMAP(integrated,n.neighbors = 5, dims = 1:20)

PCAPlot(integrated) 
integrated <- FindNeighbors(integrated, dims = 1:10) 
integrated <- FindClusters(integrated,dims = 1:10, resolution = c(0.1))#, 0.6, 0.8, 1.0, 1.4
integrated <- RunUMAP(integrated, dims = 1:5)

DimPlot(integrated, reduction = "umap") +
  scale_color_manual(values = colors)+
  plot_annotation(title = "1")

scRNA<- readRDS("ALL_seurat.rds")  
my_colors <- RColorBrewer::brewer.pal(n = 8, name = "Set2")
"#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3"
DimPlot(scRNA, 
        reduction = 'umap',                 
        group.by = 'Celltype',             
        pt.size = 1,                       

        label = T,
        label.size = 5,
        alpha = c(0.6))  +

  theme(legend.position = "none")+
  scale_color_manual(values = my_colors)

unique(scRNA$celltype)

scRNA$celltype <- as.character(scRNA$celltype)

scRNA$celltype[is.na(scRNA$celltype) | scRNA$celltype == "Smooth_muscle_cells"] <- "Fibroblast"

scRNA$celltype <- factor(scRNA$celltype)

unique(scRNA$orig.ident)

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggpubr)

sources <- list(
  source1 = "SeuratProject",
  source2 = c("GSM7475327", "GSM7475328"),
  source3 = c("GSM5645894", "GSM5645895", "GSM5645896"),
  source4 = c("GSM6112137_LC1", "GSM6112138_LC2", "GSM6112139_LC3", 
              "GSM6112140_LC4", "GSM6112141_LC5", "GSM6112142_LC6",
              "GSM6112143_LC7", "GSM6112144_LC8", "GSM6112145_LC9", 
              "GSM6112146_LC10"),
  source5 = c(
    "GSM6957590_NSCLC_PA001_sn_raw_feature_bc_matrix",
    "GSM6957591_NSCLC_PA004_sn_raw_feature_bc_matrix",
    "GSM6957592_NSCLC_PA005_sn_raw_feature_bc_matrix",
    "GSM6957593_NSCLC_PA019_sn_raw_feature_bc_matrix",
    "GSM6957594_NSCLC_PA025_sn_raw_feature_bc_matrix",
    "GSM6957595_NSCLC_PA034_sn_raw_feature_bc_matrix",
    "GSM6957596_NSCLC_PA042_sn_raw_feature_bc_matrix",
    "GSM6957597_NSCLC_PA043_sn_raw_feature_bc_matrix",
    "GSM6957598_NSCLC_PA048_sn_raw_feature_bc_matrix",
    "GSM6957599_NSCLC_PA054_sn_raw_feature_bc_matrix",
    "GSM6957600_NSCLC_PA056_sn_raw_feature_bc_matrix",
    "GSM6957601_NSCLC_PA060_sn_raw_feature_bc_matrix",
    "GSM6957602_NSCLC_PA067_sn_raw_feature_bc_matrix",
    "GSM6957603_NSCLC_PA068_sn_raw_feature_bc_matrix",
    "GSM6957604_NSCLC_PA070_sn_raw_feature_bc_matrix",
    "GSM6957605_NSCLC_PA072_sn_raw_feature_bc_matrix",
    "GSM6957606_NSCLC_PA076_sn_raw_feature_bc_matrix",
    "GSM6957607_NSCLC_PA080_sn_raw_feature_bc_matrix",
    "GSM6957608_NSCLC_PA104_sn_raw_feature_bc_matrix",
    "GSM6957609_NSCLC_PA125_sn_raw_feature_bc_matrix",
    "GSM6957610_NSCLC_PA141_sn_raw_feature_bc_matrix",
    "GSM6957613_NSCLC_STK_1_sn_raw_feature_bc_matrix",
    "GSM6957614_NSCLC_STK_3_sn_raw_feature_bc_matrix",
    "GSM6957615_NSCLC_STK_5dot1_sn_raw_feature_bc_matrix",
    "GSM6957616_NSCLC_STK_5dot2_sn_raw_feature_bc_matrix",
    "GSM6957625_NSCLC_KRAS_6_sn_raw_feature_bc_matrix",
    "GSM6957626_NSCLC_KRAS_7_sn_raw_feature_bc_matrix",
    "GSM6957627_NSCLC_KRAS_8_sn_raw_feature_bc_matrix",
    "GSM6957628_NSCLC_KRAS_17_sn_raw_feature_bc_matrix"
  ),
  source6 = c(
    "OMIX007088-01", "OMIX007088-03", "OMIX007088-07", "OMIX007088-11",
    "OMIX007088-13", "OMIX007088-14", "OMIX007088-15", "OMIX007088-18",
    "OMIX007088-25", "OMIX007088-26", "OMIX007088-27", "OMIX007088-29",
    "OMIX007088-30", "OMIX007088-32", "OMIX007088-33", "OMIX007088-34",
    "OMIX007088-35", "OMIX007088-37", "OMIX007088-39", "OMIX007088-41",
    "OMIX007088-43", "OMIX007088-51", "OMIX007088-52", "OMIX007088-55",
    "OMIX007088-56", "OMIX007088-57", "OMIX007088-58", "OMIX007088-59",
    "OMIX007088-60", "OMIX007088-61", "OMIX007088-62", "OMIX007088-64",
    "OMIX007088-67", "OMIX007088-75", "OMIX007088-76", "OMIX007088-78",
    "OMIX007088-79", "OMIX007088-80", "OMIX007088-81", "OMIX007088-82",
    "OMIX007088-85"
  )
)

scRNA$source <- NA_character_

for (i in seq_along(sources)) {
  sample_names <- sources[[i]]
  scRNA$source[scRNA$orig.ident %in% sample_names] <- names(sources)[i]
}

unassigned <- scRNA@meta.data %>% 
  filter(is.na(source)) %>% 
  distinct(orig.ident)

cell_prop <- scRNA@meta.data %>%
  group_by(source, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(source) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_stacked <- ggplot(cell_prop, aes(x = source, y = prop, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Cell Type Composition by Source",
       x = "Data Source",
       y = "Percentage",
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggplot2)
library(scales)
library(dplyr)

p_stacked <- ggplot(cell_prop, aes(x = source, y = prop, fill = celltype)) +
 
  geom_col(width = 0.65, position = position_stack(reverse = TRUE)) +
  
  geom_text(aes(label = ifelse(prop > 0.05, percent(prop, accuracy = 1), "")), 
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 3.5, color = "black",  show.legend = FALSE) +

  scale_y_continuous(expand = c(0, 0), 
                     labels = percent_format(accuracy = 1),
                     breaks = seq(0, 1, 0.2),
                     limits = c(0, 1.05)) +  
  scale_x_discrete(expand = expansion(add = 0.5)) +  

  scale_fill_manual(values = custom_colors) +
 
  theme_minimal(base_size = 12) +
  theme(
    text = element_text( color = "black"),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 15, margin = margin(r = 10), vjust = 2),
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,
      vjust = 1,
      size = 15,
      color = "black"
    ),
    axis.text.y = element_text(size = 15, color = "black"),
   
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.text = element_text(size = 15),
    legend.background = element_rect(fill = "white", color = NA),
    
    plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"), 
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  labs(
    title = "Cell Type Composition by Data Source",
    y = "Percentage of Cells"
  )

print(p_stacked)

