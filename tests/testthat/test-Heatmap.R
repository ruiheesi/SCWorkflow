CRObject <- getparamhm("TEC")

test_that("Produce heatmap and return filtered dataframe - TEC data", {
  
      heatplot <- do.call(Heatmap,CRObject)
      dev.off()
      print(heatplot$plot)
      expected.elements = c("plot","data")
      expect_setequal(names(heatplot), expected.elements)
})      

test_that("Heatmap scaled vs unscaled", {

  CRObject_test <- CRObject$object
  CRObject_test$scale.data <- FALSE
  heatplot2 <- do.call(Heatmap,CRObject)
  dev.off()
  print(heatplot$plot)
  
  #compare scaled (a) vs. nonscaled data (b) to be different
  a <- rowMeans(as.data.frame.matrix(heatplot$data)[,-1]) 
  b <- rowMeans(as.data.frame.matrix(heatplot2$data)[,-1])
  
  expect_false(isTRUE(all.equal(a, b)))

})      

test_that("Heatmap run with bad gene name", {
  
  CRObject_test <- CRObject
  CRObject_test$transcripts <- c("badgene",CRObject_test$transcripts) #Spike in with bad gene name
  expect_warning(do.call(Heatmap,CRObject_test),"^There are")
  
  #Actual warning:
  #There are 1 gene(s) absent from dataset:'badgene'. 
  #Possible reasons are that gene is not official gene symbol 
  #or gene is not highly expressed and has been filtered.
  
}) 

test_that("Heatmap run with duplicate gene names", {
  
  CRObject_test <- CRObject
  CRObject_test$transcripts <- c("Map7",CRObject_test$transcripts) #Spike in with bad gene name
  expect_warning(do.call(Heatmap,CRObject_test),"The following duplicate genes were removed: Map7")
  
}) 

test_that("Heatmap run with error for missing all genes", {
  
  CRObject_test <- CRObject
  CRObject_test$transcripts <- c("APOE","ARG1","CD38","CD3D","CD3E","CD3G")
  expect_error(do.call(Heatmap,CRObject_test),
               "No genes listed are found in dataset.")
  
})

test_that("Produce heatmap and return filtered dataframe - Chariou data", {
  
  CRObject <- getparamhm("Chariou")
  heatplot <- do.call(Heatmap,CRObject)
  dev.off()
  print(heatplot$plot)
  expected.elements = c("plot","data")
  expect_setequal(names(heatplot), expected.elements)
})    