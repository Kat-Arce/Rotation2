install.packages ("Seurat")
install.packages ("future")
install.packages ("ggplot2")
install.packages ("sctransform")             
install.packages("arrow")

library(Seurat)
library(future)
library(ggplot2)
library(sctransform)
library(arrow)

data.dir <- "/path/to/main/Xenium/data/folder"

#Seurat needs a csv.gz to read the transcript file but the current 10X Xenium output only has a parquet file. So the parquet needs to be converted prior to loading the Xenium object.
transcripts <- read_parquet("path/transcripts.parquet")
write.csv(transcripts, gzfile(file.path(data.dir, "transcripts.csv.gz")), row.names = FALSE)
file.exists(file.path(data.dir, "transcripts.csv.gz"))

#Load the Xenium data
path <- "/path/to/main/Xenium/data/folder"
xenium.obj <- LoadXenium(path, fov = "fov")

# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

#QC plot
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)

#To run on my mac I increased max size for globals to 8 GB (or more if needed)
options(future.globals.maxSize = 8 * 1024^3)  # Set to 8GB

#Normalize the data. Can do SCTransform or Normalize+Scale Data, whichever memory allows
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium") 

xenium.obj <- NormalizeData(xenium.obj) 
xenium.obj <- ScaleData(xenium.obj, assay = "Xenium")

#Below generates clusters needed for UMAPs + Feature Plots
xenium.obj <- RunPCA(xenium.obj, npcs = 30, features = rownames(xenium.obj))
xenium.obj <- RunUMAP(xenium.obj, dims = 1:30)
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:30)
xenium.obj <- FindClusters(xenium.obj, resolution = 0.3)

#Save xenium obj
saveRDS(xenium.obj, file = "/path/where/to/save/xenium_object.rds")

#How to load xenium obj after saving
xenium.obj <- readRDS("/path/where/to/save/xenium_object.rds")

#Below runs the view of the UMAP
DimPlot(xenium.obj) + ggtitle("UMAP Plot of X")

#UMAP Feature Plots for transcripts of interest, replace with your fave
FeaturePlot(xenium.obj, features = c("Wasf1")) 

#Image Feature Plots for transcripts of interest, replace with your fave
ImageFeaturePlot(xenium.obj, features = c("Wasf1"), cols = c("white", "red"))

#used this function to show how many cells per cluster there are to explain weird graphs for certain clusters
table(Idents(xenium.obj))

###END###
