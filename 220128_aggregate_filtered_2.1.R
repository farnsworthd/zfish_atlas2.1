###january 28, 2022 -- to be run on talapas..
library(Seurat)

###define path to matrices and sample names
#/projects/acmillerlab/drf/runs/atlas2.0/v4.3.2_millerLab_v3/outs/
#/projects/acmillerlab/drf/R/atlas2.0/

Atlas_2 <- file.path("/projects/acmillerlab/drf/runs/atlas2.0/v4.3.2_millerLab_v3/outs")
dataset_list <- c("1a","1b","2a","2b","3a","3b","4a","4b","5a","5b","5c","6a","6b","7a","7b")
merge_list <- c("1a","1b","2a","2b","3a","3b","4a","4b","5a","5b","5c","6a","6b","7a","7b")

##create seurat objects from each sample matrix
for (dataset in dataset_list){
  print(dataset)
  data <- Read10X(data.dir = file.path(Atlas_2, dataset, "outs/filtered_feature_bc_matrix"))
  atlas2.1 <- CreateSeuratObject(counts = data, project = paste0("sample",dataset), min.cells = 3, min.features = 200)
  assign(paste0("sample",dataset), atlas2.1)
}

##merge seurat objects to unify atlas object and save as .rds
#atlas2.1 <- merge(sample7a, y = c(sample7b), add.cell.ids = dataset_list, project = "atlas2.1")
atlas2.1 <- merge(sample1a, y = c(sample1b, sample2a, sample2b, sample3a, sample3b, sample4a, sample4b, sample5a, sample5b, sample5c, sample6a, sample6b, sample7a, sample7b), add.cell.ids = dataset_list, project = "atlas2.1")
print("merged")

##### seurat analysis of merged atlas object
atlas2.1[["percent.mt"]] <- PercentageFeatureSet(atlas2.1, pattern = "^mt-")
atlas2.1 <- subset(atlas2.1, subset = nFeature_RNA > 200 & nCount_RNA < 100000 & percent.mt < 25)

# store mitochondrial percentage in object meta data
atlas2.1 <- PercentageFeatureSet(atlas2.1, pattern = "^mt-", col.name = "percent.mt")

#run sctransform
atlas2.1 <- SCTransform(atlas2.1, vars.to.regress = "percent.mt", verbose = T)
print("sct_complete")
atlas2.1 <- RunPCA(atlas2.1, features = VariableFeatures(object = atlas2.1), npcs = 150)
atlas2.1 <- FindNeighbors(atlas2.1, dims = 1:115)
atlas2.1 <- FindClusters(atlas2.1, resolution = c(24.0))
atlas2.1 <- RunUMAP(atlas2.1, dims = 1:115)
saveRDS(atlas2.1, file = "/projects/acmillerlab/drf/R/atlas2.0/rds/atlas2.1_115PC_res24.rds")
print("clustering complete")

atlas2.1@misc$markers <- FindAllMarkers(atlas2.1, only.pos = TRUE, min.pct = 0.3)
colnames(atlas2.1@misc$markers)[2]  <- "avg_logFC"
saveRDS(atlas2.1, file = "/projects/acmillerlab/drf/R/atlas2.0/rds/atlas2.1_115PC_res24_markers.rds")
print("all markers complete")


