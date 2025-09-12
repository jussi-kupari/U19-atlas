library(Seurat)
library(ggplot2)

seu <- readRDS("Seurat.integr.rds")

# Trained with full mouse integrated atlas from Krauter et al.
pollock_predictions_u19 <- read.csv("U19_predictions_full_model.csv")

# Rename first column to "cell"
names(pollock_predictions_u19)[1] <- "cell"

# Select columns 
start_col <- which(names(pollock_predictions_u19) == "predicted_cell_type_probability")
end_col <- which(names(pollock_predictions_u19) == "probability.TRPM8.2")
pollock_predictions <- pollock_predictions_u19[, c(1, start_col:end_col)]

# Convert rownames to column
meta <- seu@meta.data
meta[["cell"]] <- rownames(meta)
meta <- meta[, c("cell", "Berlin.3")]
names(meta)[2] <- "celltype"

# Left join
meta <- merge(meta, pollock_predictions, by = "cell", all.x = TRUE)
meta <- meta[, c("celltype", names(meta)[5:length(names(meta))])]

names(meta) <- sub("\\.$", "", names(meta))
names(meta) <- sub("probability\\.", "", names(meta))

# Add row ID for reshape
meta$id <- seq_len(nrow(meta))

# Reshape to long format
cols_to_pivot <- names(meta)[!names(meta) %in% c("celltype", "id")]

meta_long <- reshape(
  meta,
  direction = "long",
  varying = cols_to_pivot,
  v.names = "score",
  timevar = "name",
  times = cols_to_pivot,
  idvar = c("celltype", "id")
)

# Plot violins for cell type prediction scores
prediction_violins <-
  meta_long |>
  ggplot(aes(celltype, score, fill = celltype)) +
  geom_violin(scale = "width", alpha = 0.4) +
  geom_jitter(size = 0.1, alpha = 0.05) +
  geom_boxplot(alpha = 1,
               width = 0.1,
               outlier.shape = NA) +
  facet_wrap( ~ name,
              strip.position = "right",
              ncol = 1,
              scales = "free_y") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey99", colour = "grey20"),
    strip.text.y.right = element_text(angle = 0),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )




# Reshape data to wide matrix
mean_scores <- aggregate(score ~ celltype + name, data = meta_long, FUN = mean)
mean_scores$temp_id <- seq_len(nrow(mean_scores))
wide_data <- reshape(mean_scores,
                     direction = "wide",
                     idvar = "celltype",
                     timevar = "name",
                     v.names = "score")

names(wide_data) <- gsub("^score\\.", "", names(wide_data))
wide_data$temp_id <- NULL

rownames(wide_data) <- wide_data$celltype
wide_data$celltype <- NULL
pollock_matrix <- as.matrix(wide_data)

# Order rows to approx match mouse cell types
rows_order <- c(
  "Ab-LTMR.LGI2",
  "Ab-LTMR.ETV1",
  "Ab-LTMR.NSG2",
  "Ab-LTMR.CCKAR",
  "ATF3",
  "C-LTMR.CDH9",
  "C-NP.MRGPRX1/GFRA2" ,
  "C-NP.MRGPRX1/MRGPRX4",
  "C-NP.SST",
  "A-LTMR.TAC3",
  "C-PEP.ADORA2B",
  "C-PEP.TAC1/CHRNA3",
  "C-PEP.TAC1/CACNG5",
  "A-PEP.CHRNA7/SLC18A3",
  "A-PEP.NTRK3/S100A16",
  "A-PEP.KIT",
  "A-PEP.SCGN/ADRA2C",
  "A-Propr.HAPLN4",
  "A-Propr.EPHA3",
  "A-Propr.PCDH8",
  "C-Thermo.RXFP1",
  "C-Thermo.TRPM8"
) 

pollock_mat <- pollock_matrix[rows_order, ]

# Plot heatmap of cell type predictions with z-scored rows 
prediction_heatmap <-
  pheatmap::pheatmap(
    pollock_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "column",
    number_color = "white",
    fontsize_number = 8,
    color = colorRampPalette(c("blue", "white", "red"))(100)
  )