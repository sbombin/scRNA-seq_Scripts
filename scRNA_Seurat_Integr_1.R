library(Seurat)
library(readr)
library(dplyr)
library(Matrix)
library(future)
#library(clustree)

### the order of samples in metadata file should match the order in the filepath file

ReadObject <- function(name){ readRDS(paste0(data_path, name, ".RDS")) }
SaveObject <- function(object, name){ saveRDS(object, paste0(data_path, name, ".RDS")) }

data_path <- "seurat_results/"
fig_path <- "seurat_results/"

SaveFigure <- function(plots, name, type = "pdf", width, height, res){
  if(type == "pdf") {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

#### Find automatically an optimal Number PCAs to explain variation
MIN_CUM_SD <- 90 ## Percent of variance to be captured by PCs.
## Get number of PCs required to capture xx% of variability                 
fgetNPCs <- function(obj,MIN_CUM_SD=90){
  sds <- Stdev(obj,reduction="pca")
  cumsds <- 0
  for(i in 1:length(sds)){
    cumsds <- cumsds + sds[i]
    if(cumsds >= MIN_CUM_SD){
      return(i)
    }
  }
}



metadata <- read_csv("metadata.csv")

filepath <- read_csv("filepath.txt", col_names = FALSE)
filepath <-  filepath$X1

data <- lapply(filepath, Read10X)

id <- metadata$id

# name samples in the list
names(data) <- id

# create empty list
object.list <- list()

# add seurat objects to list
for (x in seq_along(data)) {
  object.list[[x]] <- CreateSeuratObject(counts = data[[x]], project = id[x], min.cells = 3, min.features = 200) 
}

# name samples in the list
names(object.list) <- id

## MT ratio
for (x in seq_along(object.list)) {
  object.list[[x]][["percent.mt"]] = PercentageFeatureSet(object.list[[x]], pattern = "^MT-")
}

## find upper filtering threshold
score.count <- lapply(object.list, function(x) mean(x$nCount_RNA) + 3*sd(x$nCount_RNA))
score.feature <- lapply(object.list, function(x) mean(x$nFeature_RNA) + 3*sd(x$nFeature_RNA))
# other option can be median absolute deviation (MAD)
#mad.count <- lapply(object.list, function(i) 10^(median(log10(i$nCount_RNA)) + 2.5*mad(log10(i$nCount_RNA))))
#mad.feature <- lapply(object.list, function(i) 10^(median(log10(i$nFeature_RNA)) + 2.5*mad(log10(i$nFeature_RNA))))

## Collect summary (optional)
nCount <- lapply(object.list, function(x) summary(x$nCount_RNA))
nFeature <- lapply(object.list, function(x) summary(x$nFeature_RNA))
nMT <- lapply(object.list, function(x) summary(x$percent.mt))
nCount
nFeature 
nMT

## Plot QC
lapply(seq_along(object.list), function(i) {
  p1 <- VlnPlot(object.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  SaveFigure(p1, paste0('QC_1_',names(object.list)[i]), width = 10, height = 10)})

lapply(seq_along(object.list), function(i) {
  p1 <- FeatureScatter(object.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  SaveFigure(p1, paste0('QC_2_',names(object.list)[i]), width = 10, height = 6)
  p2 <- FeatureScatter(object.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  SaveFigure(p2, paste0('QC_3_',names(object.list)[i]), width = 10, height = 6)})

#Add categories  
for (x in seq_along(object.list)) {
  object.list[[x]][["treatment"]] <- metadata$treatment[x]
  object.list[[x]][["timepoint"]] <- metadata$timepoint[x]
  object.list[[x]][["study"]] <- metadata$study[x]
}

## save object before filtering
SaveObject(object.list, "object_list_raw")

# Filter; standart percent.mt < 5 but it can filter out too many cells
for (x in seq_along(object.list)) {
  object.list[[x]] <- subset(object.list[[x]], subset = nFeature_RNA > 200 & nFeature_RNA < score.feature[[x]] & percent.mt < 15 &
                               nCount_RNA >500 & nCount_RNA < score.count[[x]])
}

## enable Parallelization
plan("multicore", workers = 28) ## if in Rstudio plan("multisession")
plan()
options(future.globals.maxSize = 100 * 1024 ^ 3)

## Run normalization and find features separately for each sample
object.list <- lapply(X = object.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

## select features that are repeatedly variable across datasets for integration
# integrate samples to correct for batch effect
features <- SelectIntegrationFeatures(object.list = object.list)
object.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)
object.combined <- IntegrateData(anchorset = object.anchors)


## Run the standard workflow for visualization and clustering
object.combined <- ScaleData(object.combined, verbose = FALSE)
object.combined <- RunPCA(object.combined, npcs = 50, verbose = FALSE)

## Find optimal number of PCs
object.combined <- JackStraw(object.combined, num.replicate = 100, dims = 50)
object.combined <- ScoreJackStraw(object.combined, dims = 1:50)
p <- JackStrawPlot(object.combined, dims = 1:50)
SaveFigure(p, "JackStrawPlot_PC50", width = 10, height = 8)
p <- ElbowPlot(object.combined, ndims = 50)
SaveFigure(p, "ElbowPlot_PC50", width = 10, height = 8)
npcs <- fgetNPCs(object.combined,90)
npcs


SaveObject(object.combined, "object_integrated")



