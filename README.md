# NormComapre

# Comparing 3 different Normalization stratagies on samples with variable seq depth (78-462M reads per Sample)


```{r}
library(Seurat)
library(patchwork)
#Load Seurat object named TCs
# Perform SCTransform normalization
TCs_sct <- SCTransform(TCs, verbose = TRUE) #TCs is the Seurat object I have for my samples
normalized_depth_sct <- colSums(GetAssayData(TCs_sct, slot = "data", assay = "SCT"))
TCs_sct <- AddMetaData(TCs_sct, metadata = normalized_depth_sct, col.name = "normalized_depth_sct")

# Perform LogNormalize normalization
TCs_log <- NormalizeData(TCs, normalization.method = "LogNormalize", scale.factor = median(TCs$nCount_RNA))
normalized_depth_log <- colSums(GetAssayData(TCs_log, slot = "data", assay = "RNA"))
TCs_log <- AddMetaData(TCs_log, metadata = normalized_depth_log, col.name = "normalized_depth_log")

# Perform CLR normalization
TCs_clr <- NormalizeData(TCs, normalization.method = "CLR", scale.factor = median(TCs$nCount_RNA))
normalized_depth_clr <- colSums(GetAssayData(TCs_clr, slot = "data", assay = "RNA"))
TCs_clr <- AddMetaData(TCs_clr, metadata = normalized_depth_clr, col.name = "normalized_depth_clr")

# Prepare data for plotting
plot_data_sct <- data.frame(Samples = TCs_sct@meta.data$orig.ident, NormalizedDepth = TCs_sct@meta.data$normalized_depth_sct)
plot_data_log <- data.frame(Samples = TCs_log@meta.data$orig.ident, NormalizedDepth = TCs_log@meta.data$normalized_depth_log)
plot_data_clr <- data.frame(Samples = TCs_clr@meta.data$orig.ident, NormalizedDepth = TCs_clr@meta.data$normalized_depth_clr)

# Generate individual boxplots
plot_sct <- ggplot(plot_data_sct, aes(x = Samples, y = NormalizedDepth)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = format(..y.., digits = 2, nsmall = 2)), vjust = -0.5, color = 'cyan3', size = 3.5) +
  theme_minimal() +
  labs(title = "Normalized Sequencing Depth: SCTransform", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_log <- ggplot(plot_data_log, aes(x = Samples, y = NormalizedDepth)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = format(..y.., digits = 2, nsmall = 2)), vjust = -0.5, color = 'cyan3', size = 3.5) +
  theme_minimal() +
  labs(title = "Normalized Sequencing Depth: LogNormalize", x = "", y = "Normalized Sequencing Depth (UMI Count per Cell)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_clr <- ggplot(plot_data_clr, aes(x = Samples, y = NormalizedDepth)) +
  geom_boxplot() +
  stat_summary(fun = median, geom = "text", aes(label = format(..y.., digits = 2, nsmall = 2)), vjust = -0.5, color = 'cyan3', size = 3.5) +
  theme_minimal() +
  labs(title = "Normalized Sequencing Depth: CLR", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine plots using patchwork
combined_plot <- plot_sct + plot_log + plot_clr + plot_layout(ncol = 1)

# Display the combined plot
print(combined_plot)
```
# This will generate a combined boxplot

![nc](https://github.com/user-attachments/assets/3ec03dfa-dc16-4364-86d0-68534fdf5ffa)

