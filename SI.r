library(Seurat)


# Load the ChP datasets

tumors <- c("C1", "C3", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12")
tumor_types <- vector(mode = "list", length = 14)
tumor_types <- c("Normal", "Normal", "CPP", "CPP", "CPP", "aCPP", "CPC", "CPC", "CPC", "CPC", "CPC", "CPP", "CPP", "CPC")
names(tumor_types) <- tumors
status_list <- vector(mode = "list", length = 14)
status_list <- c('Normal', 'Normal', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor','Tumor', 'Tumor')
names(status_list) <- tumors
tumors <- c("C1", "T3","T7")
#tumors <- c( "C3", "T1", "T3", "T5", "T7", "T11")
data.folder = "E:/Tony/Matrices/"
choroid.list <- vector("list", length(tumors))
i <- 1
for (tumor in tumors){
  file_loc <- paste0(data.folder,tumor,"/raw_feature_bc_matrix")
  data <- Read10X(data.dir = file_loc)
  # Initialize the Seurat object with the raw (non-normalized data).
  object <- CreateSeuratObject(counts = data, project = "BatchCorrect", min.cells = 3, min.features = 200)
  #filter cells
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
  object
  object$batch <- tumor
  object$type <- tumor_types[tumor]
  object$status <- status_list[tumor]
  choroid.list[i] <- object
  print(tumor)
  i <- i+1
}

#Normalize data
for (i in 1:length(choroid.list)) {
  choroid.list[[i]] <- NormalizeData(choroid.list[[i]], verbose = TRUE)
  choroid.list[[i]] <- FindVariableFeatures(choroid.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = TRUE)
}

#Integrate samples
choroid.anchors <- FindIntegrationAnchors(object.list = choroid.list, dims = 1:30)
choroid.integrated <- IntegrateData(anchorset = choroid.anchors, dims = 1:30)

DefaultAssay(choroid.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
choroid.integrated <- ScaleData(choroid.integrated, verbose = FALSE)
choroid.integrated <- RunPCA(choroid.integrated, npcs = 30, verbose = FALSE)
choroid.integrated <- RunUMAP(choroid.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(choroid.integrated, reduction = "umap", group.by = "batch")
p2 <- DimPlot(choroid.integrated, reduction = "umap", group.by = "ident", label = TRUE, 
              repel = TRUE) + NoLegend()
plot_grid(p1, p2)

