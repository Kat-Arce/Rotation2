library(Seurat)
library(future)
library(ggplot2)
library(sctransform)

#All_markers DE table creates comparison against each cluster vs all clusters
#p_val : p-value (un-adjusted)
#avg_log2FC : log fold-change of the average expression between the two groups. 
#     Positive values indicate that the feature is more highly expressed in the first group.

#pct.1 : The percentage of cells where the feature is detected in the first group
#pct.2 : The percentage of cells where the feature is detected in the second group (all other clusters)
#p_val_adj : Adjusted p-value, based on Bonferroni correction using all features in the dataset.

xenium.obj <- readRDS("/path/to/xenium_object.rds")

all_markers <- FindAllMarkers(xenium.obj)
View(all_markers)

write.csv(all_markers, "/path/to/save/all_markers.csv", row.names = FALSE)

#How to re-access
all_markers <- read.csv("/path/to/all_markers.csv")

#Lognormalized gene counts for all genes across all cells in each cluster with Seurat pseudo bulking
total_exp <- AggregateExpression(
  xenium.obj,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1,
  verbose = TRUE,)

total_exp_df <- as.data.frame(as.matrix(total_exp$Xenium))
View(total_exp_df)


# Step 1: Add a column with the row sums
total_exp_df$RowSum <- rowSums(total_exp_df)

# Step 2: Filter for specific transcripts
genes_to_filter <- c("Zap70", "Itga1", "Cd4", "Tnf", "Cd3e", 
                     "Tuba1a", "Tubb5", "Cd40", "Lck", "Lat", 
                     "Mapk1", "Cd19", "Itga7", "Itga4", "Itgb5", 
                     "Itgax", "Itgal", "Itgam", "Actn1", "Actn4", 
                     "Actr3", "Actr2", "Tubb5", "Tubg1", "Il2", "Itgad",
                     "Cxcl10", "Cd8a")

# Filter the data to include only these genes
filtered_data <- total_exp_df[rownames(total_exp_df) %in% genes_to_filter, ]
View(filtered_data)

#Make it a csv
write.csv(filtered_data, "/pat/to/save/aggregated_expression_filtered.csv", row.names = TRUE)
