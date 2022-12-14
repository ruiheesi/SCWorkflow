test_that("Load testing dataset", {
  # Seurat_Object <- readRDS('./otherData/PostFilterQC.rds')
  # Seurat_Object <- readRDS('./fixtures/SO_moduleScore.rds')
  
  Seurat_Object <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  
  
  
  PCA_and_Normalization_out <- PCA_and_Normalization(Seurat_Object,
                                                     vars_to_regress = c('percent.mt'),
                                                     vars_to_plot = c('percent.mt','nCount_RNA'),
                                                     npcs = 30,
                                                     nfeatures = 2000,
                                                     low_cut = 1,
                                                     high_cut = 8,
                                                     low_cut_disp = 1,
                                                     high_cut_disp = 100000,
                                                     selection_method = "vst",
                                                     doJackStraw = FALSE,
                                                     JackStraw_dims = 5,
                                                     methods_PCA = c("Elbow","Marchenko-Pastur"),
                                                     var_threshold = 0.1,
                                                     imageType = "png"
                                 )

  
  
  
  # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_PCAnorm.pdf"),
         # PCA_and_Normalization_out$plot)
  # plot(PCA_and_Normalization_out$plot)
  
  
  
  ### Test parameters
  
  # create output
  expected.elements = c("so","plot")
  expect_setequal(names(PCA_and_Normalization_out), expected.elements)
  # SO contains object same length as input
  expect_equal(length(PCA_and_Normalization_out$so),length(Seurat_Object)) 
  # figure slot is a grob
  expect_equal(class(PCA_and_Normalization_out$plot)[3], 'grob')
  # SO slot contains data
  expect( object.size(PCA_and_Normalization_out$so[[1]]@assays$RNA@counts),'> 0' )
  # plot slot contains data
  expect( object.size(PCA_and_Normalization_out$plot),'> 0' )
  
})


# library(devtools)
# document()
# load_all()
# test_active_file()

# vars_to_regress = c('percent.mt')
# vars_to_plot = c('percent.mt','nCount_RNA')
# npcs = 30
# nfeatures = 2000
# low_cut = 1
# high_cut = 8
# low_cut_disp = 1
# high_cut_disp = 100000
# selection_method = 'vst'
# doJackStraw = FALSE
# JackStraw_dims = 5
# methods_PCA = c('Elbow','Marchenko-Pastur')
# var_threshold = 0.1
# imageType = 'png'

# saveRDS(PCA_and_Normalization_out$so,'/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/PCAnorm.rds')

# ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_PCAnormgrob.pdf"), grobs)
