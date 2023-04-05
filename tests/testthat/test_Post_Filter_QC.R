
for (data in c('TEC','Chariou','pbmc-single','NSCLC-Multi')) {	
  
  test_that(paste0("Test Post Filter and QC - Standard (",data," dataset)"), {
    
    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds"))) 
    post.filter.QC.out <- postFilterQC(object$object,	
                                       image.type = 'png'
    )
    
    
    
    # create output
    expected.elements = c("so","plot")
    expect_setequal(names(post.filter.QC.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(post.filter.QC.out$so),length(object$so))
    # figure slot is a grob
    expect_equal(class(post.filter.QC.out$plot)[3], 'grob')
    # SO slot contains data
    expect( object.size(post.filter.QC.out$so[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(post.filter.QC.out$plot),'> 0' )
    
  })
  
}
