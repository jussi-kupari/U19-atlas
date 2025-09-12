## -----------------------------------------------------------------------------
library(Seurat)


## -----------------------------------------------------------------------------
mouse <- readRDS("integrated_atlas_humanized.rds") 


## -----------------------------------------------------------------------------
human <- readRDS("Seurat.integr.rds")


## -----------------------------------------------------------------------------
mouse@meta.data$Species <- "Mouse"
mouse@meta.data$CellType <- paste(mouse@meta.data$Species,
                                  mouse@meta.data$cell_type,
                                  sep = "-")

human@meta.data$Species <- "Human"
human@meta.data$CellType <- paste(human@meta.data$Species,
                                  human@meta.data$Berlin.3,
                                  sep = "-")

## -----------------------------------------------------------------------------
merged <- 
  merge(human, mouse) |> 
  JoinLayers()


## -----------------------------------------------------------------------------
merged@meta.data |> 
  count(orig.ident) |> 
  arrange(n)


## -----------------------------------------------------------------------------
# Merge some orig.idents to get group sizes to > 30

# Create a lookup vector for mappings
lookup <- c(
  "P26402_1004"    = "P26402_1003_4",
  "P26402_1003"    = "P26402_1003_4",
  "hDRGwc2"        = "hDRGwc2_p3",
  "hDRGp3"         = "hDRGwc2_p3",
  "hDRG2_Chat-cre" = "hDRG_combo1",
  "hDRG1_Chat-cre" = "hDRG_combo1",
  "hDRG_092917_2"  = "hDRG_combo1",
  "hDRG_092917_1"  = "hDRG_combo1",
  "m2r1"           = "m1_2r1",
  "m1r1"           = "m1_2r1"
)

# Apply lookup
merged@meta.data$data_set <- lookup[merged@meta.data$orig.ident]
merged@meta.data$data_set[is.na(merged@meta.data$data_set)] <-
  merged@meta.data$orig.ident[is.na(merged@meta.data$data_set)]

# Handle the SeuratProject case by using the Study column
seuratproject_mask <- merged@meta.data$orig.ident == "SeuratProject"

merged@meta.data$data_set[seuratproject_mask] <-
  merged@meta.data$Study[seuratproject_mask]


## -----------------------------------------------------------------------------
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$data_set)
merged


## -----------------------------------------------------------------------------
# run standard anlaysis workflow
merged <-
  NormalizeData(merged) |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA()


## -----------------------------------------------------------------------------
merged <-
  FindNeighbors(merged, dims = 1:30, reduction = "pca") |>
  FindClusters(resolution = 2, cluster.name = "unintegrated_clusters")


## -----------------------------------------------------------------------------
merged <- RunUMAP(merged, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


## -----------------------------------------------------------------------------
merged <-
  IntegrateLayers(
    object = merged,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca", 
    k.weight = 35,
    verbose = TRUE
  )

# re-join layers after integration
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

merged <- FindNeighbors(merged, reduction = "integrated.cca", dims = 1:30)
merged <- FindClusters(merged, resolution = 1)


## -----------------------------------------------------------------------------
merged <- RunUMAP(merged, dims = 1:30, reduction = "integrated.cca")


## -----------------------------------------------------------------------------
saveRDS(merged, file = "seurat_mouse_human_drg_integrated_cca.rds")
