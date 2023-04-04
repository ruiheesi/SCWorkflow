
for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {

  test_that(paste0("Test Post Filter and QC - Standard (",data," dataset)"), {
         
  Seurat_Object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
  Post_filter_QC_out <- Post_filter_QC(object$so,
                                 imageype = 'png'
                                 )
                            
  
  
  # create output
  expected.elements = c("so","plot")
  expect_setequal(names(Post_filter_QC_out), expected.elements)
  # SO contains object same length as input
  expect_equal(length(Post_filter_QC_out$so),length(Seurat_Object$so))
  # figure slot is a grob
  expect_equal(class(Post_filter_QC_out$plot)[3], 'grob')
  # SO slot contains data
  expect( object.size(Post_filter_QC_out$so[[1]]@assays$RNA@counts),'> 0' )
  # plot slot contains data
  expect( object.size(Post_filter_QC_out$plot),'> 0' )
  
})

}

