## --------------------------------------------------------------------------------------------------------
library(Seurat)


## --------------------------------------------------------------------------------------------------------
u19 <- readRDS("Seurat.integr.rds")


## --------------------------------------------------------------------------------------------------------
Idents(u19) <- "Berlin.3"
u19 <- BuildClusterTree(u19, dims = 1:50, reduction = "pca")


## --------------------------------------------------------------------------------------------------------
tree <- Tool(u19, slot = "BuildClusterTree")
unrooted_tree <- ape::unroot(tree)


## --------------------------------------------------------------------------------------------------------
plot(
  unrooted_tree,
  type = "unrooted",
  lab4ut = "axial",
  edge.width = 2,
  cex = 1,
  no.margin = TRUE
)


## --------------------------------------------------------------------------------------------------------
tree_distances <-
  ape::cophenetic.phylo(u19@tools[["BuildClusterTree"]])

clust_distances <- hclust(dist(tree_distances))
tree_dist_matrix <- as.matrix(tree_distances)


## --------------------------------------------------------------------------------------------------------
pheatmap::pheatmap(
  tree_dist_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  treeheight_row = 0,
  treeheight_col = 0,
  color = viridis::viridis(100)
)  

