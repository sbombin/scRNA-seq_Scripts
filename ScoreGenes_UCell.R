source('HelpFUNs.R')
library(UCell)

fig_path <- "UCell_out"

object <- readRDS("seuratobj.RDS")
DefaultAssay(object) <- "RNA"

## find files with pathway genes and load as list
file_list <- list.files( pattern = "*.grp")
file_list 
df_list <- lapply(file_list,
                  FUN = function(files) {
                    read.table(files, header=TRUE, quote="\"")
                  })

## extract gene set (pathway) names
names <- lapply(df_list, function(x) names(x))
names <- unlist(names)
## extract genes for each pathway and add to list
path <- lapply(df_list, unlist, use.names=FALSE)
names(path) <- names

## change gene names to upper case if needed
#path <- lapply(path, toupper)
## extract all gene names present in seurat object
all.genes <- rownames(object)

## remove pathway genes that are NOT in seurat dataset
path <- lapply(path, function(x) x[x %in% all.genes])
## or remove pathway genes that are present in seurat dataset
#path <- lapply(path, function(x) x[!x %in% all.genes])

## run UCell scoring
ucell <- AddModuleScore_UCell(object, features = path, maxRank = 2000, slot = "data")
signature.names <- paste0(names(path), "_UCell")

## create and save  Feature plots and Violin plots
lapply(seq_along(signature.names), function(i) {
  p1 <- VlnPlot(ucell, features = signature.names[i]) + theme(legend.position = 'none')
  SaveFigure(p1, paste0('VlnPlot_',signature.names[i]), width = 8, height = 8)})

lapply(seq_along(signature.names), function(i) {
  p1 <- VlnPlot(ucell, features = signature.names[i], split.by = "group") 
  SaveFigure(p1, paste0('VlnPlot2_',signature.names[i]), width = 10, height = 8)})

lapply(seq_along(signature.names), function(i) {
  p1 <- FeaturePlot(ucell, features = signature.names[i]) 
  SaveFigure(p1, paste0('FeaturePlot_',signature.names[i]), width = 8, height = 8)})

