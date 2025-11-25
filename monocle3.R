
cds <- as.cell_data_set(scRNA)
cds
## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)
## Add gene names into CDS
#cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(pbmc3[["RNA"]])
# to get cell metadata
colData(cds)
# to gene metdata
fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)

# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have

# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 
Idents(pbmc3)<-pbmc3$RNA_snn_res.0.2
list_cluster <- pbmc3@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- pbmc3@reductions$umap@cell.embeddings

# plot
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'RNA_snn_res.0.2',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.before.trajectory# | cluster.names

cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- align_cds(cds)

# marker finding
marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'RNA_snn_res.0.2',   #RNA_snn_res.0.3
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)


# ...4. Order the cells in pseudotime -------------------

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 4]))#, root_cells = colnames(cds[,clusters(cds) == 4])
?order_cells
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(RNA_snn_res.0.2, monocle3_pseudotime, median), fill = RNA_snn_res.0.2)) +
  geom_boxplot()

# seurat_clusters
# ...5. Finding genes that change as a function of pseudotime --------------------

deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)
print(cds)

deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

pr_deg_ids <- row.names(subset(deg_bcells, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$RNA_snn_res.0.2)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

Track_genes_sig<-c("SERPINB3", "FMO2","SERPINB4" ,"KRT5" ,"VMO1" , "BPIFB1", "MSMB" , "GABRP" , "RARRES1", "CAPN13",  
                   "CXCL6" ,"BPIFA1" ,"SCGB1A1" , "FBLN1" , "LYNX1" ,"SAA2")

Track_genes_sig<-c( "MTDH" ,"KLF4" ,"CCL20", "HSF1" ,"TMPRSS4", "IGFBP7","EFEMP1","BAG3","CCND1","PFN2","EGFL7","ZFAS1", "GSN", "BOP1", "CXCL14" ,
                    "SMC3" ,"TIMP3" ,"NR2F2" ,"S100P" ,"LGR4" ,"ELK3" )

Track_genes_sig<-c("SOX2" , "IGF1R" , "CLDN1" ,"SERPINF1" ,"MUC4" ,"AQP5" ,"SLC9A3R1","ZFP36" ,"ST6GAL1" )

#基因表达趋势图
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="RNA_snn_res.0.2", #RNA_snn_res.0.3
                         min_expr=0.5, ncol = 2)
#FeaturePlot图
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=TRUE,
           label_cell_groups=TRUE,  label_leaves=TRUE)
