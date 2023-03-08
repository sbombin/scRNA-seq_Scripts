### Convenience functions

## save plots
SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 300)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

## save RDS file
SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

## import RDS file
ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}

### Find automatically Optimal Number PCAs to explain variation
## Percent of variance to be captured by PCs.
MIN_CUM_SD <- 90 
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


### save nested list in different folders
## borrowed from MollahLab/gen3DNet

write_all_to_disk <- function(name, object, folder) {
  if (class(object) == "data.frame") {
    final_path = file.path(folder, paste(name,".csv",sep=""))
    write.csv(
      object,
      final_path
    )
  }
  else if (class(object) == "list") {
    new_path <- file.path(folder, name)
    dir.create(new_path, showWarnings=FALSE)
    lapply(
      1:length(object),
      function(idx) {
        object_name <- names(object)[idx]
        if (length(object_name) > 0) {
          write_all_to_disk(
            object_name,
            object[[idx]],
            new_path
          )
        }
      }
    )
  }
}
