

for (data in c('TEC','Chariou','PBMC_Single','NSCLC_Multi')) {
  
  test_that(
    paste0("Test Combine & Renormalize - Standard (",data," dataset)"), {
      
      
      data.run <- getParamCN(data)
      
      combine.renormalize.out <- do.call(combineNormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plots")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
      # SO slot contains data
      expect( nrow(combine.renormalize.out$object@assays$RNA@counts),'> 0' )
      # plot slot contains data
      expect( object.size(combine.renormalize.out$plot),'> 0' )
      
    })
}




for (data in c('TEC')) {
  
  test_that(
    paste0("Test Combine & Renormalize - only.var.genes (",data," dataset)"), {
      
      
      data.run <- getParamCN(data)
      
      data.run$only.var.genes = TRUE
      combine.renormalize.out <- do.call(combineNormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plots")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
      # SO slot contains data
      expect( nrow(combine.renormalize.out$object@assays$RNA@counts),'> 0' )
      # plot slot contains data
      expect( object.size(combine.renormalize.out$plot),'> 0' )
      
    })
} 



for (data in c('TEC')) {
  
  test_that(
    paste0("Test Combine & Renormalize - SCT level - sample (",data," dataset)"), {
      
      
      data.run <- getParamCN(data)
      
      data.run$SCT.level = 'Sample'
      data.run$object = data.run$object
      combine.renormalize.out <- do.call(combineNormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plots")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
      # SO slot contains data
      expect( nrow(combine.renormalize.out$object@assays$RNA@counts),'> 0' )
      # plot slot contains data
      expect( object.size(combine.renormalize.out$plot),'> 0' )
      
    })
} 


for (data in c('TEC')) {
  
  test_that(
    paste0("Test Combine & Renormalize - exclude.sample (",data," dataset)"), {
      
      data.run <- getParamCN(data)
      # data.run$input=data.run$input[c(1,2)]
      data.run$exclude.sample=data.run$object[1]%>%names
      combine.renormalize.out <- do.call(combineNormalize, data.run)
      
      
      # create output
      expected.elements = c("object","plots")
      expect_setequal(names(combine.renormalize.out), expected.elements)
      # figure slot is a grob
      expect_equal(class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
      # SO slot contains data
      expect( nrow(combine.renormalize.out$object@assays$RNA@counts),'> 0' )
      # plot slot contains data
      expect( object.size(combine.renormalize.out$plot),'> 0' )
      # sample is removed
      expect_false(
        data.run$exclude.sample%in%
          unique(combine.renormalize.out$object@meta.data$orig.ident)
                  )
      
    })
}


for (data in c('TEC')) {
  
  test_that(
    paste0("Test Combine & Renormalize - selection.method = mean.var.plot (",
           data," dataset)"), {
             
             
             data.run <- getParamCN(data)
             data.run$selection.method = "mean.var.plot"
             combine.renormalize.out <- do.call(combineNormalize, data.run)
             
             
             # create output
             expected.elements = c("object","plots")
             expect_setequal(names(combine.renormalize.out), expected.elements)
             # figure slot is a grob
             expect_equal(class(combine.renormalize.out$plots$TSNE)[3], 'ggplot')
             # SO slot contains data
             expect( nrow(combine.renormalize.out$object@assays$RNA@counts),'> 0' )
             # plot slot contains data
             expect( object.size(combine.renormalize.out$plot),'> 0' )
             
           })
}






