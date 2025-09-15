## --------------------------------------------------------------------------------------------------------
library(Seurat)


## --------------------------------------------------------------------------------------------------------
merged <- readRDS("seurat_mouse_human_drg_integrated_cca.rds") 


## --------------------------------------------------------------------------------------------------------
# Split UMAP by species
DimPlot(merged, reduction = "umap", group.by = c("Species"), shuffle = TRUE, split.by = "Species")


## --------------------------------------------------------------------------------------------------------
# Split UMAP by orig.ident
DimPlot(
  merged,
  reduction = "umap",
  group.by = c("orig.ident"),
  split.by = "Species",
  shuffle = TRUE
) 


## --------------------------------------------------------------------------------------------------------
# Split UMAP by Study

# First combine Human and Mouse study identities
merged@meta.data <-
  transform(merged@meta.data, 
            Study = ifelse(Species == "Human", Study, dataset))
DimPlot(
  merged,
  reduction = "umap",
  group.by = c("Study"),
  split.by = "Species",
  shuffle = TRUE
) 


## --------------------------------------------------------------------------------------------------------
# Split UMAP by CellType
DimPlot(
  merged,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  group.by = c("CellType"),
  split.by = "Species"
) + NoLegend()


## --------------------------------------------------------------------------------------------------------
# Build cluster tree
Idents(merged) <- "CellType"
merged <- BuildClusterTree(merged, dims = 1:30, reduction = "integrated.cca")


## --------------------------------------------------------------------------------------------------------
# Plot tree facing right
PlotClusterTree(merged, direction = "rightwards")


## --------------------------------------------------------------------------------------------------------
# Get tree and unroot it
tree <- Tool(merged, slot = "BuildClusterTree")
unrooted_tree <- ape::unroot(tree)


## --------------------------------------------------------------------------------------------------------
# Plot the rootless dendrogram
 plot(
    unrooted_tree,
    type = "unrooted",
    lab4ut = "axial",
    edge.width = 2,
    cex = 1,
    no.margin = TRUE
  )


## --------------------------------------------------------------------------------------------------------
# Create a matrix of distances from the tree
tree_distances_merged <-
  ape::cophenetic.phylo(merged@tools[["BuildClusterTree"]])

clust_distances_merged <- hclust(dist(tree_distances_merged))
tree_dist_df_merged <- as.data.frame(tree_distances_merged)
tree_dist_matrix_merged <- as.matrix(tree_dist_df_merged)


## --------------------------------------------------------------------------------------------------------
# Plot a heatmap from the distance matrix
pheatmap::pheatmap(
  tree_dist_matrix_merged,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  treeheight_row = 0,
  treeheight_col = 0,
  color = viridis::viridis(100)
)  
