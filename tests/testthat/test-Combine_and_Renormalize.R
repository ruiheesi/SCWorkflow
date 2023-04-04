

for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  
  test_that(
    paste0("Test Combine & Renormalize - Standard (",data," dataset)"), {
      
      
      data.run <- getParamCR(data)
      combine.renormalize.out <- do.call(combineRenormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plot")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plot)[3], 'grob')
      # SO slot contains data
      expect( nrow(combine.renormalize.out$object@assays$RNA@counts),'> 0' )
      # plot slot contains data
      expect( object.size(combine.renormalize.out$plot),'> 0' )
      
    })
}



for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  
  test_that(
    paste0("Test Combine & Renormalize - NoRegress (",data," dataset)"), {
      
      
      data.run <- getParamCR(data)
      data.run$vars.to.regress = NULL
      combine.renormalize.out <- do.call(combineRenormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plot")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plot)[3], 'grob')
      # SO slot contains data
      expect( nrow(combine.renormalize.out$object@assays$RNA@counts),'> 0' )
      # plot slot contains data
      expect( object.size(combine.renormalize.out$plot),'> 0' )
      
    })
}


for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  
  test_that(
    paste0("Test Combine & Renormalize - only.var.genes (",data," dataset)"), {
      
      
      data.run <- getParamCR(data)
      data.run$only.var.genes = TRUE
      combine.renormalize.out <- do.call(combineRenormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plot")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plot)[3], 'grob')
      # SO slot contains data
      expect( nrow(combine.renormalize.out$object@assays$RNA@counts),'> 0' )
      # plot slot contains data
      expect( object.size(combine.renormalize.out$plot),'> 0' )
      
    })
} 



for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  
  test_that(
    paste0("Test Combine & Renormalize - exclude.sample (",data," dataset)"), {
      
      data.run <- getParamCR(data)
      data.run$exclude.sample=object[1]%>%names
      combine.renormalize.out <- do.call(combineRenormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plot")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plot)[3], 'grob')
      # SO slot contains data
      expect( nrow(combine.renormalize.out$object@assays$RNA@counts),'> 0' )
      # plot slot contains data
      expect( object.size(combine.renormalize.out$plot),'> 0' )
      
    })
}


for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  
  test_that(
    paste0("Test Combine & Renormalize - selection.method = mean.var.plot (",
           data," dataset)"), {
             
             
             data.run <- getParamCR(data)
             data.run$selection.method = "mean.var.plot"
             combine.renormalize.out <- do.call(combineRenormalize, data.run)
             
             
             # create output
             expected.elements = c("object","plot")
             expect_setequal(names(combine.renormalize.out), expected.elements)
             # figure slot is a grob
             expect_equal(class(combine.renormalize.out$plot)[3], 'grob')
             # SO slot contains data
             expect( nrow(combine.renormalize.out$object@assays$RNA@counts),'> 0' )
             # plot slot contains data
             expect( object.size(combine.renormalize.out$plot),'> 0' )
             
           })
}



for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  
  test_that(
    paste0("Test Combine & Renormalize - integrate.data (",data," dataset)"), {
      
      
      data.run <- getParamCR(data)
      data.run$integrate.data = TRUE
      combine.renormalize.out <- do.call(combineRenormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plot")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plot)[3], 'grob')
      # SO slot contains data
      expect( nrow(combine.renormalize.out$object@assays$RNA@counts),'> 0' )
      # plot slot contains data
      expect( object.size(combine.renormalize.out$plot),'> 0' )
      
    })
}


