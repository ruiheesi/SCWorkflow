for (data in c('TEC','Chariou','NSCLC_Multi')) {#,'PBMC_Single')) {
  
  test_that(paste0("Test Filter and QC - Standard (",data," dataset)"), {
    
    data.run <- getParamRaw(data)
    Raw.out <- do.call(processRawData, data.run)
    
    
    # create output
    expected.elements = c("object","plots")
    expect_setequal(names(Raw.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(Raw.out$object),length(data.run$input))
    # figure slot is a ggplot
    expect_equal(class(Raw.out$plots[[1]])[2], 'ggplot')
    # SO slot contains data
    expect( object.size(Raw.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(Raw.out$plots),'= 0' )
    
  })
  
}


################################################################

for (data in c('TEC')) {
  
  test_that(paste0("Test Filter and QC - Regex (",data," dataset)"), {
    
    data.run <- getParamRaw(data)
    data.run$keep = TRUE
    data.run$file.filter.regex = "1|2"
    Raw.out <- do.call(processRawData, data.run)
    
    
    # create output
    expected.elements = c("object","plots")
    expect_setequal(names(Raw.out), expected.elements)
    # SO contains object same length as input
    expect_false(equals(length(Raw.out$object),length(data.run$input)))
    # figure slot is a ggplot
    expect_equal(class(Raw.out$plots[[1]])[2], 'ggplot')
    # SO slot contains data
    expect( object.size(Raw.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(Raw.out$plots),'>0' )
    
  })
  
}


