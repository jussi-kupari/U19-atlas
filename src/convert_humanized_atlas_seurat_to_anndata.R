library(Seurat)
library(SeuratDisk)

seu <- readRDS("krauter_integrated_mouse_drg_atlas_humanized.rds")

# Convert to V3/4/Assay structure
seu <- scCustomize::Convert_Assay(seu, convert_to = "V3")

DefaultAssay(seu) <- "RNA"

seu@assays$integrated <-
  seu@assays$RNA$scale.data <-
  seu@assays$SCT <-
  seu@graphs$SCT_nn <-
  seu@graphs$SCT_snn <-
  seu@misc$scPred <- NULL

h5Seu <- paste("humanized_mouse_integrated_atlas", "h5Seurat", sep = ".")
SaveH5Seurat(seu, filename = h5Seu)
Convert(h5Seu, dest = "h5ad", assay = "RNA")
