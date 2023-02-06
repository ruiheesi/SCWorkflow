downsample_SO=function(IN,OUT){
  library(Seurat)
  
SO=readRDS(IN)
downsampled.obj <- SO[sample(rownames(SO), size =3000, replace=F), sample(colnames(SO), size =2000, replace=F)]

# dir.create(paste0(pathOut,"tests/testthat/fixtures/",x), showWarnings = FALSE)
saveRDS(downsampled.obj,file=OUT)

}



# ########################################################################################
# ## CITEseq data did not work with previous function. Found this function that works for CITEseq
# 
# # https://rdrr.io/github/bimberlabinternal/CellMembrane/src/R/SeuratUtils.R
# DownsampleSeurat <- function(seuratObj, targetCells, subsetFields = NULL, seed = GetSeed()) {
#   if (!is.null(seed)) {
#     set.seed(seed)
#   }
#   
#   if (!is.null(subsetFields)){
#     for (subsetField in subsetFields) {
#       if (!(subsetField %in% names(seuratObj@meta.data))) {
#         stop(paste0('Field not found in seuratObj: [', subsetField, ']'))
#       }
#     }
#   }
#   
#   print(paste0('Subsetting, original cells: ', ncol(seuratObj)))
#   cellsToRetain <- c()
#   if (is.null(subsetFields)) {
#     toSample <- min(targetCells, ncol(seuratObj))
#     if (toSample != targetCells) {
#       print(paste0('There are only ', ncol(seuratObj), ' cells. Will not downsample.'))
#       return(seuratObj)
#     }
#     
#     cellsToRetain <- sample(colnames(seuratObj), toSample, replace = F)
#   } else {
#     groupVals <- (seuratObj@meta.data %>% tidyr::unite("x", subsetFields, remove = FALSE))$x
#     names(groupVals) <- colnames(seuratObj)
#     counts <- table(groupVals)
#     print(paste0('Unique values: ', paste0(unique(names(counts)), collapse = ',')))
#     for (val in unique(names(counts))) {
#       availBarcodes <- colnames(seuratObj)[groupVals == val]
#       toSample <- min(targetCells, length(availBarcodes))
#       if (toSample != targetCells) {
#         print(paste0('There are only ', length(availBarcodes), ' cells available for group: ', val, '. Will retain all cells for this group.'))
#         cellsToRetain <- c(cellsToRetain, availBarcodes)
#       } else {
#         cellsToRetain <- c(cellsToRetain, sample(availBarcodes, targetCells, replace = F))
#       }
#     }
#   }
#   
#   print(paste('Total cells retained ', length(cellsToRetain), ' of ', ncol(seuratObj)))
#   seuratObj <- subset(seuratObj, cells = cellsToRetain)
#   
#   if (ncol(seuratObj) != length(cellsToRetain)) {
#     stop(paste('Incorrect number of cells retained!, was: ', ncol(seuratObj), ', expected: ', length(cellsToRetain)))
#   }
#   
#   print(paste0('Final cells: ', ncol(seuratObj)))
#   
#   return(seuratObj)
# }
# 
# 
# ############################################
# downsample_CITE_SO=function(IN,OUT){
#   
#   SO=readRDS(IN)
#   downsampled.obj <- DownsampleSeurat(SO,targetCells=2000,seed = 42)
#   # downsampled.obj=downsampled.obj[sample(rownames(downsampled.obj), size =3000, replace=F), ]
#   if('Protein'%in%names(downsampled.obj@assays)){
#     genes=c(sample(rownames(downsampled.obj), size =3000-length(rownames(downsampled.obj@assays$Protein)), replace=F),rownames(downsampled.obj@assays$Protein))%>%unique()
#   }else{
#     genes=c(sample(rownames(downsampled.obj), size =3000, replace=F))%>%unique()
#   }
#   downsampled.obj <- subset(downsampled.obj,features = genes)
#   
#   # dir.create(paste0(pathOut,"tests/testthat/fixtures/",x), showWarnings = FALSE)
#   saveRDS(downsampled.obj,file=OUT)  
# }
# 
# 
# 
# 
# 
# ########################################################################################
# ########################################################################################
# ########################################################################################
# 
# 
# for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
#   # data='NSCLC_Multi'
#   
#     if (data=='TEC') {org='Mouse'
#     }else{org='Human'}
#     
# 
#     localFilePaths=list.files(paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/",data,"/h5files"),".h5",full.names = T)#%>%as.list
#   
#       print(paste0("Filter_and_QC - ",data))
#     filter_qc_out <- Filter_and_QC(localFilePaths,
#                                    organism = org,
#                                    rename = F,
#                                    New_Sample_Names = c("Sample_1", "Sample_2"),
#                                    mincells = 10,
#                                    mingenes = 500,
#                                    complexity = 0.6,
#                                    MAD_genes_value = 3,
#                                    MAD_mitoch_value = 3,
#                                    minUMI = 500,
#                                    MAD_gene = TRUE,
#                                    MAD_mitoch = TRUE,
#                                    maxgenes = 2500,
#                                    maxmitoch = 10,
#                                    Filter_VDJ_Genes = FALSE,
#                                    Keep = TRUE,
#                                    File_Filter_Regex = "",
#                                    Split_H5 = FALSE,
#                                    Protein = FALSE,
#                                    Cell_hash = FALSE,
#                                    imageType = "png",
#                                    plot_histogram = FALSE
#     )%>%suppressMessages()%>%suppressWarnings()
#     
#     
#     downsampled.obj=lapply(filter_qc_out$so, function(x){DownsampleSeurat(x,targetCells=1000,seed = 42)})
#     if('Protein'%in%names(downsampled.obj@assays)){
#       genes=c(sample(rownames(downsampled.obj), size =1500-length(rownames(downsampled.obj@assays$Protein)), replace=F),rownames(downsampled.obj@assays$Protein))%>%unique()
#     }else{
#       genes=c(sample(rownames(downsampled.obj), size =1500, replace=F))%>%unique()
#     }
#     downsampled.obj <- subset(downsampled.obj,features = genes)    
#     filter_qc_out$so=downsampled.obj
# 
#     OUT=paste0('/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/',data,'/',data,'_Filtered_SO_downsample.rds')
#     saveRDS(filter_qc_out,file=OUT)    
#     
#     
#     print(paste0("PCA Norm - ",data))
#     PCA_and_Normalization_out <- PCA_and_Normalization(filter_qc_out$so,
#                                                        vars_to_regress = c('percent.mt'),
#                                                        vars_to_plot = c('percent.mt','nCount_RNA'),
#                                                        npcs = 30,
#                                                        nfeatures = 2000,
#                                                        low_cut = 1,
#                                                        high_cut = 8,
#                                                        low_cut_disp = 1,
#                                                        high_cut_disp = 100000,
#                                                        selection_method = "vst",
#                                                        doJackStraw = FALSE,
#                                                        JackStraw_dims = 5,
#                                                        methods_PCA = c("Elbow","Marchenko-Pastur"),
#                                                        var_threshold = 0.1,
#                                                        imageType = "png"
#     )
# 
#     saveRDS(PCA_and_Normalization_out,
#             file=paste0('/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/fixtures/',data,'/',data,'_Filtered_PCA_Norm_SO_downsample.rds'))    
#     
#     rm("PCA_and_Normalization_out","filter_qc_out",'downsampled.obj','OUT','localFilePaths')
#     
#     } 
# 
# 
# 
# ########################################################################################
# ########################################################################################
# ########################################################################################


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



# Git Token
# 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
