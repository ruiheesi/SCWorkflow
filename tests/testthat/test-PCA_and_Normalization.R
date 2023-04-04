for (data in c('TEC','Chariou','pbmc-single','nsclc-multi')) {

  test_that(
    paste0("Test PCA and Normalization - Standard (",data," dataset)"),
            {
    
    
    data.run <- getParamPN(data)
    pca.normalization.out <- do.call(pcaAndNormalization, data.run)
    
  
  ### Test parameters
  # create output
  expected.elements = c("object","plot")
  expect_setequal(names(pca.normalization.out), expected.elements)
  # SO contains object same length as input
  expect_equal(length(pca.normalization.out$object),length(object$object)) 
  # figure slot is a grob
  expect_equal(class(pca.normalization.out$plot)[3], 'grob')
  # SO slot contains data
  expect(object.size(
    pca.normalization.out$object[[1]]@assays$RNA@counts),'>0')
  # plot slot contains data
  expect( object.size(pca.normalization.out$plot),'> 0' )
  
})

}




################################################################################

for (data in c('TEC','Chariou','pbmc-single','nsclc-multi')) {

  test_that(
    paste0("Test PCA and Normalization - NoRegress (",data," dataset)"), 
            {

    
    data.run <- getParamPN(data)
    data.run$vars.to.regress = c()
    pca.normalization.out <- do.call(pcaAndNormalization, data.run)    
    
    
    ### Test parameters
    # create output
    expected.elements = c("object","plot")
    expect_setequal(names(pca.normalization.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(pca.normalization.out$object),length(object))
    # figure slot is a grob
    expect_equal(class(pca.normalization.out$plot)[3], 'grob')
    # SO slot contains data
    expect(object.size(
      pca.normalization.out$object[[1]]@assays$RNA@counts),'>0')
    # plot slot contains data
    expect(object.size(pca.normalization.out$plot),'> 0' )

  })

}



################################################################################

for (data in c('TEC','Chariou','pbmc-single','nsclc-multi')) {

  test_that(
    paste0("Test PCA and Normalization - Jackstraw-elbow (",data," dataset)"), 
    {

    data.run <- getParamPN(data)
    data.run$jackstraw = TRUE
    data.run$methods.pca = c("Elbow")
    
    pca.normalization.out <- do.call(pcaAndNormalization, data.run)   


    ### Test parameters
    # create output
    expected.elements = c("object","plot")
    expect_setequal(names(pca.normalization.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(pca.normalization.out$object),length(object))
    # figure slot is a grob
    expect_equal(class(pca.normalization.out$plot)[3], 'grob')
    # SO slot contains data
    expect(object.size(
      pca.normalization.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(pca.normalization.out$plot),'> 0' )

  })

}


################################################################################

for (data in c('TEC','Chariou','pbmc-single','nsclc-multi')) {

  test_that(
    paste0("Test PCA and Normalization-Jackstraw-Marchenko-Pastu (",
           data," dataset)"), {

    data.run <- getParamPN(data)
    data.run$jackstraw = TRUE
    data.run$methods.pca = c("Marchenko-Pastur")
    pca.normalization.out <- do.call(pcaAndNormalization, data.run)   
    

    ### Test parameters
    # create output
    expected.elements = c("object","plot")
    expect_setequal(names(pca.normalization.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(pca.normalization.out$object),length(object))
    # figure slot is a grob
    expect_equal(class(pca.normalization.out$plot)[3], 'grob')
    # SO slot contains data
    expect(object.size(
      pca.normalization.out$object[[1]]@assays$RNA@counts),'>0')
    # plot slot contains data
    expect(object.size(pca.normalization.out$plot),'> 0' )

  })

}


################################################################################

for (data in c('TEC','Chariou','pbmc-single','nsclc-multi')) {

  test_that(
    paste0("Test PCA and Normalization - selection- mean.var.plot (",
           data," dataset)"), {

    data.run <- getParamPN(data)
    data.run$selection.method = "mean.var.plot"
    pca.normalization.out <- do.call(pcaAndNormalization, data.run)   


    ### Test parameters
    # create output
    expected.elements = c("object","plot")
    expect_setequal(names(pca.normalization.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(pca.normalization.out$object),length(object))
    # figure slot is a grob
    expect_equal(class(pca.normalization.out$plot)[3], 'grob')
    # SO slot contains data
    expect(object.size(
      pca.normalization.out$object[[1]]@assays$RNA@counts),'>0')
    # plot slot contains data
    expect( object.size(pca.normalization.out$plot),'> 0' )

  })

}




################################################################################

for (data in c('TEC','Chariou','pbmc-single','nsclc-multi')) {

  test_that(
    paste0("Test PCA and Normalization - selection- dispersion (",
                   data," dataset)"), {

    data.run <- getParamPN(data)
    data.run$selection.method = "dispersion"
    pca.normalization.out <- do.call(pcaAndNormalization, data.run)   


    ### Test parameters
    # create output
    expected.elements = c("object","plot")
    expect_setequal(names(pca.normalization.out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(pca.normalization.out$object),length(object))
    # figure slot is a grob
    expect_equal(class(pca.normalization.out$plot)[3], 'grob')
    # SO slot contains data
    expect(object.size(
      pca.normalization.out$object[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(pca.normalization.out$plot),'> 0' )

  })

}


