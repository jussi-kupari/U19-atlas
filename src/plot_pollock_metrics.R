library(ggplot2)

metrics <- read.csv("full_atlas_training_metrics.csv") 


# Reshape to long format
cols_to_pivot <- names(metrics)[!names(metrics) %in% c("epoch")]

metrics_long <-
  reshape(
    metrics,
    direction = "long",
    varying = cols_to_pivot,
    v.names = "value",
    timevar = "metric",
    times = cols_to_pivot,
    idvar = c("epoch")
)

metrics_long |>
  subset(metric %in% c("train_accuracy", "val_accuracy")) |>
  ggplot(aes(epoch, value, color = metric)) +
  geom_line() +
  cowplot::theme_cowplot() +

metrics_long |>
  subset(metric %in% c("train.loss", "val.loss")) |>
  ggplot(aes(epoch, value, color = metric)) +
  geom_line() +
  cowplot::theme_cowplot()
