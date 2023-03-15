library(Seurat)
library(sigclust)
#library(MultiK)
library(ggdendro)

library(doParallel)
library(foreach)

object <- readRDS("seuratobj.Sham-BDL.scale_all.RDS")

start.time <- Sys.time()

#multik <- MultiK(object, reps=120, nPC =22, resolution = seq(0.1, 1, 0.5))
#saveRDS(multik,"multik.sham-bdl.RDS")
#p <- DiagMultiKPlot(multik$k, multik$consensus)
#ggsave('Optimal_Ks.sham-bdl.pdf', plot = p, scale = 1.5)

#clusters <- getClusters(object, 16)
#saveRDS(clusters,"clusters.sham-bdl.RDS")

RunSigClust <- function(x1, x2, l1, l2) {
  suppressPackageStartupMessages(library(sigclust))
  sig.dat <- as.matrix(t(cbind(x1, x2)))
  l1 <- rep(1, length(l1))
  l2 <- rep(2, length(l2))
  sig.label <- c(l1, l2)
  sig <- sigclust(x = sig.dat,
                  nsim = 100,
                  nrep = 1, labflag = 1, # uses known label
                  label = sig.label,
                  icovest = 2)

  return(sig)
}

CalcSigClustParal <- function(object, clusters) {
  suppressPackageStartupMessages(library(sigclust))
  
  hvg <- VariableFeatures(object=object)
  norm.hvg <- object@assays$RNA@data[hvg, ]
  ClustAssign <- as.character(clusters)
  n <- length(unique(ClustAssign))
  #pval <- matrix(NA, ncol = n, nrow = n)
  #rownames(pval) <- c(paste("cluster", c(0: (n-1)), sep = ""))
  #colnames(pval) <- c(paste("cluster", c(0: (n-1)), sep = ""))
  
  # set up data to run in SigClust
  x.list <- list()
  l.list <- list()
  for (i in 1: n) {
    x.list[[paste0("cluster", i-1)]] <- as.matrix(norm.hvg[, ClustAssign == i-1])
    l.list[[i]] <- ClustAssign[ClustAssign == i-1]
  }
  
  # run SigClust
  #pval <- data.frame()
  pval <- foreach (i= 1: (n-1), .packages='sigclust') %:%
    foreach (j=(i+1): n) %dopar% {
      RunSigClust(x1 = x.list[[i]],
                                x2 = x.list[[j]],
                                l1 = l.list[[i]],
                                l2 = l.list[[j]])@pval
    }
  pval
  
}


parallel::detectCores()
registerDoParallel(parallel::detectCores()-2) # set number of cores to available -2 


pval <- CalcSigClustParal(object, clusters=object$RNA_snn_res.0.6)
saveRDS(pval,"pval.list.RDS")

#p <- PlotSigClust(object, clusters=object$RNA_snn_res.0.6, pval)
#ggsave('Significant_Clusters.pdf', plot = p, scale = 1.8)


end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
