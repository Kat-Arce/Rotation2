library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(Seurat)

#### Read and process objects in one at a time due to device memory limitations ####
Y1xenium.obj <- readRDS("path/to/object.rds")
Y2xenium.obj <- readRDS("path/to/object.rds")
O1xenium.obj <- readRDS("path/to/object.rds")
O2xenium.obj <- readRDS("path/to/object.rds")

#Used GetAssay to figure out how to access the correct metadata. 
GetAssay(Y1xenium.obj, assay = "Xenium")


###################################################################################
#### How to Extract the scaled/normalized expression matrix with z scores per cell ####
expression_matrix <- GetAssayData(Y1xenium.obj, assay = "Xenium", layer = "scale.data")
expression_matrix <- GetAssayData(Y2xenium.obj, assay = "Xenium", layer = "scale.data")
expression_matrix <- GetAssayData(O1xenium.obj, assay = "Xenium", layer = "scale.data")
expression_matrix <- GetAssayData(O2xenium.obj, assay = "Xenium", layer = "scale.data")

# Calculate the mean z score of each gene
gene_means <- rowMeans(expression_matrix)

# Convert to a data frame for easier viewing 
gene_mean_df <- data.frame(Gene = rownames(expression_matrix), MeanExpression = gene_means)
View(gene_mean_df)
write.csv(gene_mean_df, "path/to/save.csv", row.names = FALSE)
write.csv(gene_mean_df, "path/to/save.csv", row.names = FALSE)
write.csv(gene_mean_df, "path/to/save.csv", row.names = FALSE)
write.csv(gene_mean_df, "path/to/save.csv", row.names = FALSE)

###################################################################################
#### How to Extract the log normalized expression matrix ####
expression_matrix <- GetAssayData(Y1xenium.obj, assay = "Xenium", layer = "data")
expression_matrix <- GetAssayData(Y2xenium.obj, assay = "Xenium", layer = "data")
expression_matrix <- GetAssayData(O1xenium.obj, assay = "Xenium", layer = "data")
expression_matrix <- GetAssayData(O2xenium.obj, assay = "Xenium", layer = "data")

# Calculate the mean expression of each gene
gene_means <- rowMeans(expression_matrix)

# Convert to a data frame for easier viewing 
gene_mean_df <- data.frame(Gene = rownames(expression_matrix), MeanExpression = gene_means)
View(gene_mean_df)

write.csv(gene_mean_df, "path/to/save.csv", row.names = FALSE)
write.csv(gene_mean_df, "path/to/save.csv", row.names = FALSE)
write.csv(gene_mean_df, "path/to/save.csv", row.names = FALSE)
write.csv(gene_mean_df, "path/to/save.csv", row.names = FALSE)
###################################################################################
####  How to extract counts ####
#The matrix is too large to extract each cell x transcript count 
#so I summed all cells per transcript 
raw_counts <- GetAssayData(O2xenium.obj, assay = "Xenium", layer = "counts")
gene_sums <- rowSums(raw_counts)
gene_sums_df <- data.frame(
  Gene = rownames(raw_counts),
  TotalCounts = gene_sums
)


#I added a total row to sum all of the counts within the sample
total_sum <- sum(gene_sums)
gene_sums_df <- rbind(gene_sums_df, data.frame(Gene = "Total", TotalCounts = total_sum))

View(gene_sums_df)

write.csv(gene_sums_df, "path/to/save.csv", row.names = FALSE)
write.csv(gene_sums_df, "path/to/save.csv", row.names = FALSE)
write.csv(gene_sums_df, "Dpath/to/save.csv", row.names = FALSE)
write.csv(gene_sums_df, "path/to/save.csv", row.names = FALSE)
###################################################################################

#### How to filter out genes of interest for z scores, can be applied to means but change the code ####
genes_to_filter <- c("Zap70", "Itga1", "Cd4", "Cd3e",
                     "Tuba1a", "Tubb5", "Cdc42", "Cd40",
                     "Lck", "Lat", "Mapk1", "Cd19", "Itga7", "Cd44", "Itga4", "Itgb5", 
                     "Itgax", "Itgal", "Itgam", "Actn1", 
                     "Actn4", "Actr3", "Actr2", "Tubg1","Itgad", "Cxcl10", "Cd8a" )



# Filter the data frame to include only the specified genes
filtered_genes_df <- subset(gene_mean_df, Gene %in% genes_to_filter)
View(filtered_genes_df)
write.csv(filtered_genes_df, "path/to/save.csv", row.names = FALSE)
write.csv(filtered_genes_df, "path/to/save.csv", row.names = FALSE)
write.csv(filtered_genes_df, "path/to/save.csv", row.names = FALSE)
write.csv(filtered_genes_df, "path/to/save.csv", row.names = FALSE)
###################################################################################

#### How to filter and normalize counts ####
# List of genes to filter, including "Total"
total_genes_to_filter <- c("Zap70", "Itga1", "Cd4", "Cd3e",
                     "Tuba1a", "Tubb5", "Cdc42", "Cd40",
                     "Lck", "Lat", "Mapk1", "Cd19", "Itga7", "Cd44", "Itga4", "Itgb5", 
                     "Itgax", "Itgal", "Itgam", "Actn1", 
                     "Actn4", "Actr3", "Actr2", "Tubg1", "Itgad", "Cxcl10", "Cd8a", "Total")

# Make it a data frame
filtered_gene_sums_df <- subset(gene_sums_df, Gene %in% total_genes_to_filter)

# Extract the total sum value from the "Total" row
total_count <- filtered_gene_sums_df[filtered_gene_sums_df$Gene == "Total", "TotalCounts"]
head(total_count)

# Add a column for normalized counts
filtered_gene_sums_df$NormalizedCounts <- filtered_gene_sums_df$TotalCounts / total_count
View(filtered_gene_sums_df)

write.csv(filtered_gene_sums_df, "path/to/save.csv", row.names = FALSE)
write.csv(filtered_gene_sums_df, "path/to/save.csv", row.names = FALSE)
write.csv(filtered_gene_sums_df, "path/to/save.csv", row.names = FALSE)
write.csv(filtered_gene_sums_df, "path/to/save.csv", row.names = FALSE)
###################################################################################
##END##


