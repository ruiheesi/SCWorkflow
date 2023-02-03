
downsample_SO=function(IN,OUT){
  library(Seurat)
  
SO=readRDS(IN)
downsampled.obj <- SO[sample(rownames(SO), size =3000, replace=F), sample(colnames(SO), size =2000, replace=F)]

# dir.create(paste0(pathOut,"tests/testthat/fixtures/",x), showWarnings = FALSE)
saveRDS(downsampled.obj,file=OUT)

}


## TEC
IN='/rstudio-files/ccbr-data/data/singlecell/TEC/TEC_Combine_and_Renormalize_SO.rds'
OUT='/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/TEC/TEC_Combine_and_Renormalize_SO_downsample.rds'
downsample_SO(IN,OUT)
# cp /rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/TEC/TEC_Combine_and_Renormalize_SO_downsample.rds \
# /rstudio-files/ccbr-data/data/singlecell/TEC

IN='/rstudio-files/ccbr-data/data/singlecell/TEC/TEC_Cell_Types_SingleR_SO.rds'
OUT='/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/TEC/TEC_CellTypesSingleR_SO_downsample.rds'
downsample_SO(IN,OUT)
# cp /rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/TEC/TEC_CellTypesSingleR_SO_downsample.rds
# /rstudio-files/ccbr-data/data/singlecell/TEC



## Chariou
IN='/rstudio-files/ccbr-data/data/singlecell/Chariou/Chariou_Combine_and_Renormalize_SO.rds'
OUT='/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/Chariou/Chariou_Combine_and_Renormalize_SO_downsample.rds'
downsample_SO(IN,OUT)
# cp /rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/Chariou/Chariou_Combine_and_Renormalize_SO_downsample.rds \
# /rstudio-files/ccbr-data/data/singlecell/Chariou


IN='/rstudio-files/ccbr-data/data/singlecell/Chariou/Chariou_Cell_Types_SingleR_SO.rds'
OUT='/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/Chariou/Chariou_Cell_Types_SingleR_SO_downsample.rds'
downsample_SO(IN,OUT)
# cp /rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/Chariou/Chariou_Cell_Types_SingleR_SO_downsample.rds \
# /rstudio-files/ccbr-data/data/singlecell/Chariou



## NSCLC_Single

IN='/rstudio-files/ccbr-data/data/singlecell/NSCLC_Single/NSCLCsingle_Combine_and_Renormalize_SO.rds'
OUT='/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/NSCLC_Single/NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds'
downsample_SO(IN,OUT)
# cp /rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/NSCLC_Single/NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds \
# /rstudio-files/ccbr-data/data/singlecell/NSCLC_Single


IN='/rstudio-files/ccbr-data/data/singlecell/NSCLC_Single/NSCLCsingle_Cell_Types_SingleR_SO.rds'
OUT='/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/NSCLC_Single/NSCLCsingle_Cell_Types_SingleR_SO_downsample.rds'
downsample_SO(IN,OUT)
# cp /rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/NSCLC_Single/NSCLCsingle_Cell_Types_SingleR_SO_downsample.rds \
# /rstudio-files/ccbr-data/data/singlecell/NSCLC_Single


## NSCLC_multi

IN='/rstudio-files/ccbr-data/data/singlecell/NSCLC_Multi/NSCLCmulti_Combine_and_Renormalize_SO.rds'
OUT='/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/NSCLC_Multi/NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds'
downsample_SO(IN,OUT)
# cp /rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/NSCLC_Multi/NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds \
# /rstudio-files/ccbr-data/data/singlecell/NSCLC_Multi


IN='/rstudio-files/ccbr-data/data/singlecell/NSCLC_Multi/NSCLCmulti_Cell_Types_SingleR_SO.rds'
OUT='/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/NSCLC_Multi/NSCLCmulti_Cell_Types_SingleR_SO_downsample.rds'
downsample_SO(IN,OUT)
# cp /rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/NSCLC_Multi/NSCLCmulti_Cell_Types_SingleR_SO_downsample.rds \
# /rstudio-files/ccbr-data/data/singlecell/NSCLC_Multi


## BRCA

IN='/rstudio-files/ccbr-data/data/singlecell/BRCA/BRCA_Combine_and_Renormalize_SO.rds'
OUT='/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/BRCA/BRCA_Combine_and_Renormalize_SO_downsample.rds'
downsample_SO(IN,OUT)
# cp /rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/BRCA/BRCA_Combine_and_Renormalize_SO_downsample.rds \
# /rstudio-files/ccbr-data/data/singlecell/BRCA


