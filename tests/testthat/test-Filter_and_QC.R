for (data in c('TEC','Chariou','pbmc-single','nsclc-multi')) {
  
  test_that(paste0("Test Filter and QC - Standard (",data," dataset)"), {
    
    
    data.run <- getParamFQ(data)
    filter.qc.out <- do.call(filterQC, data.run)
    
    
    # create output
    expected.elements = c("object","plot")
    expect_setequal(names(filter.qc.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(filter.qc.out$object),length(localFilePaths))
    # figure slot is a grob
    expect_equal(class(filter.qc.out$plot)[3], 'grob')
    # SO slot contains data
    expect( object.size(filter.qc.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(filter.qc.out$plot),'> 0' )
    
  })
  
}

################################################################

for (data in c('TEC','Chariou','pbmc-single','nsclc-multi')) {
  
  test_that(paste0("Test Filter and QC - plot.histogram (",data," dataset)"), {
    
    
    data.run <- getParamFQ(data)
    data.run$plot.histogram <- "TRUE"
    filter.qc.out <- do.call(filterQC, data.run)
    
    
    # create output
    expected.elements = c("object","plot")
    expect_setequal(names(filter.qc.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(filter.qc.out$object),length(localFilePaths))
    # figure slot is a grob
    expect_equal(class(filter.qc.out$plot)[3], 'grob')
    # SO slot contains data
    expect( object.size(filter.qc.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(filter.qc.out$plot),'> 0' )
    
  })
  
}


################################################################

for (data in c('TEC','Chariou','pbmc-single','nsclc-multi')) {
  # data='   } else if (data == "pbmc-single'
  
  test_that(paste0("Test Filter and QC - filter.vdj.genes (",data," dataset)"), {
    
    data.run <- getParamFQ(data)
    data.run$filter.vdj.genes <- "TRUE"
    filter.qc.out <- do.call(filterQC, data.run)
    
    
    # create output
    expected.elements = c("object","plot")
    expect_setequal(names(filter.qc.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(filter.qc.out$object),length(localFilePaths))
    # figure slot is a grob
    expect_equal(class(filter.qc.out$plot)[3], 'grob')
    # SO slot contains data
    expect( object.size(filter.qc.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(filter.qc.out$plot),'> 0' )
    
  })
  
}

################################################################

for (data in c('TEC','Chariou','pbmc-single','nsclc-multi')) {
  
  test_that(paste0("Test Filter and QC - Regex (",data," dataset)"), {
    
    data.run <- getParamFQ(data)
    data.run$keep = TRUE
    data.run$file.filter.regex = "Sample_[1-2]"
    filter.qc.out <- do.call(filterQC, data.run)
    
    
    # create output
    expected.elements = c("object","plot")
    expect_setequal(names(filter.qc.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(filter.qc.out$object),length(localFilePaths))
    # figure slot is a grob
    expect_equal(class(filter.qc.out$plot)[3], 'grob')
    # SO slot contains data
    expect( object.size(filter.qc.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(filter.qc.out$plot),'> 0' )
    
  })
  
}


