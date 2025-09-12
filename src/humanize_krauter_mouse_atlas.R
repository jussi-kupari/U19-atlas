library(Seurat)

seu <- readRDS("krauter_integrated_mouse_drg_atlas.rds")
metadata <- seu@meta.data

seu <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "counts"))

# Biomart (https://www.ensembl.org/index.html -> BioMart)
# Archive site (bottom right corner) -> choose Ensembl 102
# Choose database: Ensembl genes 102
# Mouse Genes (GRCm38.p6)
# Homologs: ->
# Gene stable ID
# Gene name

# Orthologs: [F-J] ->
# Human gene stable ID
# Human gene name
# Human homology type
# Human Gene-order conservation score
# Human orthology confidence [0 low, 1 high]
# Results as tsv -> mart_export.txt

orthos <-
  read.delim("mart_export.txt", sep = "\t") |>
  subset(
    Human.homology.type == "ortholog_one2one" &
      Human.Gene.order.conservation.score >= 90 &
      Gene.name %in% rownames(seu)
  ) |>
  subset(select = c(Gene.name, Human.gene.name))

seu <- seu[orthos$Gene.name, ]

# Order 'genes' df by filtered_matrix rownames
orthos <- orthos[match(rownames(seu), orthos$Gene.name), ]

# Change rownames to human
rownames(seu) <- make.unique(orthos$Human.gene.name)

# Create Seurat and process
seu <- CreateSeuratObject(seu)
seu@meta.data <- metadata

seu <-
  NormalizeData(seu) |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA(npcs = 100) |>
  FindNeighbors(dims = 1:20) |>
  FindClusters() |>
  RunUMAP(dims = 1:20)

seu <-
  seu |>
  harmony::RunHarmony(
    "orig.ident",
    max_iter = 50,
    reduction.use = "pca",
    reduction.save = "harmony_pca"
  )

seu <-
  FindNeighbors(seu,
                dims = seq_along(seu@reductions$harmony_pca),
                reduction = "harmony_pca") |>
  FindClusters(resolution = 0.5, verbose = FALSE) |>
  RunUMAP(
    reduction = "harmony_pca",
    dims = seq_along(seu@reductions$harmony_pca),
    verbose = FALSE,
    reduction.name = "harmony_UMAP",
    reduction.key = "harmonyUMAP_"
  )

# Save
saveRDS(seu, file = "krauter_integrated_mouse_drg_atlas_humanized.rds")
