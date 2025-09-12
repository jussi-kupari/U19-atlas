library(tidyverse)
library(Seurat)
library(SeuratDisk)

seu <- readRDS("Seurat.integr.rds")

DefaultAssay(seu) <- "RNA"

seu@assays$integrated <- NULL

seu@assays$RNA$scale.data <-
  seu@assays$SCT <-
  seu@graphs$SCT_nn <-
  seu@graphs$SCT_snn <-
  seu@misc$scPred <- NULL

# Convert to V3/4/Assay structure
seu <- scCustomize::Convert_Assay(seu, convert_to = "V3")

h5Seu <- paste("U19_atlas", "h5Seurat", sep = ".")
SaveH5Seurat(seu, filename = h5Seu)
Convert(h5Seu, dest = "h5ad", assay = "RNA")
