
downsample_SO=function(IN,OUT){
  library(Seurat)
  
SO=readRDS(IN)
downsampled.obj <- SO[sample(rownames(SO), size =3000, replace=F), sample(colnames(SO), size =2000, replace=F)]

# dir.create(paste0(pathOut,"tests/testthat/fixtures/",x), showWarnings = FALSE)
saveRDS(downsampled.obj,file=OUT)

}



IN='/rstudio-files/ccbr-data/data/singlecell/TEC/TEC_Cell_Types_SingleR_SO.rds'
OUT='/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/TEC/TEC_CellTypesSingleR_SO_downsample.rds'
downsample_SO(IN,OUT)





