test_that("Load testing dataset", {
  # Seurat_Object <- readRDS('/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/FilterQC.rds')
  # Seurat_Object <- readRDS("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/CombNorm.rds")
         
  Seurat_Object <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  
  
  Post_filter_QC_out <- Post_filter_QC(Seurat_Object,
                                 Parallelize_Computation = F, 
                                 Image_type = 'png'
                                 )
                            
  
  # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_PostFilterQC.pdf"),Post_filter_QC_out$plot)
  # plot(Post_filter_QC_out$plot)
  
  # create output
  expected.elements = c("so","plot")
  expect_setequal(names(Post_filter_QC_out), expected.elements)
  # SO contains object same length as input
  expect_equal(length(Post_filter_QC_out$so),length(Seurat_Object))
  # figure slot is a grob
  expect_equal(class(Post_filter_QC_out$plot)[3], 'grob')
  # SO slot contains data
  expect( object.size(Post_filter_QC_out$so[[1]]@assays$RNA@counts),'> 0' )
  # plot slot contains data
  expect( object.size(Post_filter_QC_out$plot),'> 0' )
  
})

# 
# library(devtools)
# document()
# load_all()
# test_active_file()


# saveRDS(Post_filter_QC_out$so,'/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/PostFilterQC.rds')
