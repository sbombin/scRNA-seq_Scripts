## load functions
source('scRNA_part1.R')
source('Save_NestList.R')

## Set directories and load seurat object
setwd("/Users/sbombin/Desktop/EICC/Projects/Frank_scRNA/")
data_path <- "/Users/sbombin/Desktop/EICC/Projects/Frank_scRNA/"
fig_path <- "/Users/sbombin/Desktop/EICC/Projects/Frank_scRNA/Plots/"

## Set number of PCs to use
pcs = 50
## Load data
object <- ReadObject("object_integrated")
DefaultAssay(object) <- "integrated"
object <- FindNeighbors(object, reduction = "pca", dims = 1:pcs)
## Check clustering with different resolutions
clustree <- FindClusters(object, resolution = c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4))
p <- clustree(clustree)
SaveFigure(p, "Clustree", width = 8, height = 10)

## optionally check the number of cells in clusters
#table(clustree$integrated_snn_res.0.4)

rm(clustree)

# Select optimal resolution 
res = 0.4
# choose grouping variable
group = "integrated_snn_res.0.4"

#make dimensional plots
object <- FindClusters(object, resolution = res)
object<- RunUMAP(object, reduction = "pca", dims = 1:pcs)
object <- RunTSNE(object, dims = 1:pcs)

p <- DimPlot(object, reduction = "umap", label = TRUE, repel = TRUE, group.by = group)
SaveFigure(p, paste0('UMAP_Clusters_',group), width = 10, height = 8)
p <- DimPlot(object, reduction = "tsne", label = TRUE, repel = TRUE, group.by = group)
SaveFigure(p, paste0('tSNE_Clusters_',group), width = 10, height = 8)

# save processed seurat object
SaveObject(object, "object_processed")

object <- ReadObject("object_processed")
### Find variable markers (genes) with DE
# Change assay
DefaultAssay(object) <- "RNA"
## optionaly create metadata column with added 'Cluster' to the cluster numbers
object$cluster <- paste0("Cluster",object$integrated_snn_res.0.4)
## Clusteers identity can be switched if needed
Idents(object) <- cluster


# function to find cluster specific markers; test.use wilcox or MAST are solid options
AllMarkers <- function(object,test){
  
  markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = test)
  # select top 10 markers for each cluster
  markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10.markers
  # save results as csv
  write.csv(top10.markers,paste0('Top10_',test,group,'.csv'), row.names = FALSE)
  write.csv(mast,paste0('AllMarkers_',test,group,'.csv'), row.names = FALSE)
}

# call function
AllMarkers(object=object,test='MAST')

## Plot markers
# create list of markers; case sensitive
markers <- list(mk1=c('FAP','COL1A2','DCN', 'COL3A1','COL6A1','PDPN'),
                mk2=c('Sox9','Sox4'))

lapply(seq_along(markers), function(i) {
  p1 <- DotPlot(object,features = markers[[i]], group.by = group) +  RotatedAxis()
  SaveFigure(p1, paste0("DotPlot_",names(markers[i]),group), width = 10, height = 8)
  p2 <- FeaturePlot(object,features = markers[[i]], reduction = "tsne")
  SaveFigure(p2, paste0("FeaturePlot_",names(markers[i]),group), width = 10, height = 8)
  p3 <- VlnPlot(object,features = markers[[i]], group.by = group) +  RotatedAxis()
  SaveFigure(p3, paste0("VlnPlot_",names(markers[i]),group), width = 10, height = 8)})
  
  
### DE groups comparison
# create vectors (or list) with control and experimental groups
# ident.1 is experimental group and ident.2 is control group; change group names according the sample names
ident.1 <- c('treat1','treat2','treat3')
ident.2 <- c('control1','control2','control3')
# combine comparison groups for naming
comparisons <- paste0(ident.1,"_vs_",ident.2)

## check number of cells in each group
ncc <- table(paste(object$cluster,object$sample, sep = "_")) %>% as.data.frame()
write.csv(ncc,'Ncell_perCluster-Sample.csv', row.names = FALSE)

## Run DE for each comparison in each cluster
# select clusters from metadata
cluster <- unique(object$cluster) %>% sort()
de.list <- list()
for (x in seq_along(comparisons)) {
  
  de.list[[x]] <- lapply(cluster, function(i) {
    FindMarkers(object, ident.1 = ident.1[x], ident.2 = ident.2[x], group.by = 'sample', subset.ident = i,
                logfc.threshold = 0.1,test.use = "wilcox",min.pct = 0.1,min.cells.group = 1) %>%
      cbind(cluster = i, .)})
}
names(de.list) <- comparisons
de.list <- lapply(de.list,setNames,cluster)

#optional function to save all comparisons in different folders
write_all_to_disk("de_results", de.list, ".")

## Make Volcano plots for DE results
library(EnhancedVolcano)

for (z in seq_along(comparisons)) {
  lapply(seq_along(cluster), function(i) {
    EnhancedVolcano(de.list[[z]][[i]], lab = row.names(de.list[[z]][[i]]), x = 'avg_log2FC', y = 'p_val_adj',
                    title = names(de.list[z]), subtitle = names(de.list[[z]][i]),
                    legendLabels=c('Not sig.','Log (base 2) FC','q-value', 'q-value & Log (base 2) FC'),
                    drawConnectors = TRUE, pCutoff = 0.01, FCcutoff = 0.5,
                    pointSize = 3.0, axisLabSize = 12, legendLabSize = 10, legendIconSize = 3.0, labSize = 3.0)
    ggsave(paste0('de_results/',names(de.list[z]),'/','Volcano_',names(de.list[[z]][i]),'.png'),dpi=300, scale =1.5)})
}

