library(UCell)

data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/"

object <- ReadObject("seuratobj.IEC.scale_all")
DefaultAssay(object) <- "RNA"
object <- RenameIdents(object, `0` = "EPL1", `1` = "MEP1", `2` = "MEP2",`3` = "TA", `4` = "MED", `5` = "IEP", `6` = "MEP3", `7` = "EPE1", `8` = "ISC",
                       `9` = "EPE2", `10` = "EPL2", `11` = "Goblet", `12` = "Tuft", `13` = "EEC", `14` = "MEP4", `15` = "Paneth")

data_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/UCell/"
fig_path <- "/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/UCell/"
setwd("/Users/sbombin/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/UCell/")

#read.table("~/Desktop/analysis_mm39.comb/Seurat_Analysis/IEC/UCell/CGAS_TARGET_GENES.v2022.1.Mm.grp", header=TRUE, quote="\"")

file_list <- list.files( pattern = "*.grp")
file_list 
df_list <- lapply(file_list,
                  FUN = function(files) {
                    read.table(files, header=TRUE, quote="\"")
                  })


names <- lapply(df_list, function(x) names(x))
names <- unlist(names)

path <- lapply(df_list, unlist, use.names=FALSE)
names(path) <- names

#path <- lapply(path, toupper)

all.genes <- rownames(object)

#a<- lapply(path, function(x) x[x %in% all.genes])
#b <- lapply(path, function(x) x[!x %in% all.genes])
#c <- as.data.frame(all.genes)

path <- lapply(path, function(x) x[x %in% all.genes])

ucell <- AddModuleScore_UCell(object, features = path, maxRank = 2000, slot = "data")
signature.names <- paste0(names(path), "_UCell")

#p <- FeaturePlot(ucell, reduction = "umap", features = signature.names, ncol = 5, order = T)
#p

lapply(seq_along(signature.names), function(i) {
  p1 <- VlnPlot(ucell, features = signature.names[i]) + theme(legend.position = 'none')
  SaveFigure(p1, paste0('VlnPlot_',signature.names[i]), width = 8, height = 8)})

lapply(seq_along(signature.names), function(i) {
  p1 <- VlnPlot(ucell, features = signature.names[i], split.by = "group") 
  SaveFigure(p1, paste0('VlnPlot2_',signature.names[i]), width = 10, height = 8)})

lapply(seq_along(signature.names), function(i) {
  p1 <- FeaturePlot(ucell, features = signature.names[i]) 
  SaveFigure(p1, paste0('FeaturePlot_',signature.names[i]), width = 8, height = 8)})

