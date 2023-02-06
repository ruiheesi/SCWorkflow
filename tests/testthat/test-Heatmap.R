CRObject <- getparamhm("TEC")

test_that("Produce heatmap and return filtered dataframe - TEC data", {
  
      output <- do.call(heatmapSC,CRObject)
      
      ggsave("output/TEC_heatmap.png",output$plot, width = 10, height = 10)
      expect_snapshot_file("output","TEC_heatmap.png")
      
      expect_type(output,"list")
      expected.elements = c("plot","data")
      expect_setequal(names(output), expected.elements)
})      

test_that("Heatmap scaled vs unscaled", {

  CRObject_test <- CRObject$object
  CRObject_test$scale.data <- FALSE
  output2 <- do.call(heatmapSC,CRObject)
  
  ggsave("output/TEC_heatmap_unscaled.png",output2$plot, width = 10, height = 10)
  expect_snapshot_file("output","TEC_heatmap_unscaled.png")
  
  #compare scaled (a) vs. nonscaled data (b) to be different
  a <- rowMeans(as.data.frame.matrix(output$data)[,-1]) 
  b <- rowMeans(as.data.frame.matrix(output2$data)[,-1])
  
  expect_false(isTRUE(all.equal(a, b)))

})      

test_that("Heatmap run with bad gene name", {
  
  CRObject_test <- CRObject
  CRObject_test$transcripts <- c("badgene",CRObject_test$transcripts) #Spike in with bad gene name
  expect_warning(do.call(heatmapSC,CRObject_test),"^There are")
  
  #Actual warning:
  #There are 1 gene(s) absent from dataset:'badgene'. 
  #Possible reasons are that gene is not official gene symbol 
  #or gene is not highly expressed and has been filtered.
  
}) 

test_that("Heatmap run with duplicate gene names", {
  
  CRObject_test <- CRObject
  CRObject_test$transcripts <- c("Map7",CRObject_test$transcripts) #Spike in with bad gene name
  expect_warning(do.call(heatmapSC,CRObject_test),"The following duplicate genes were removed: Map7")
  
}) 

test_that("Heatmap run with error for missing all genes", {
  
  CRObject_test <- CRObject
  CRObject_test$transcripts <- c("APOE","ARG1","CD38","CD3D","CD3E","CD3G")
  expect_error(do.call(heatmapSC,CRObject_test),
               "No genes listed are found in dataset.")
  
})

test_that("Produce heatmap and return filtered dataframe - Chariou data", {
  
  CRObject <- getparamhm("Chariou")
  output <- do.call(heatmapSC,CRObject)
  
  ggsave("output/Chariou_heatmap.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","Chariou_heatmap.png")
  expect_type(output,"list")
  
  expected.elements = c("plot","data")
  expect_setequal(names(output), expected.elements)
})    

test_that("Produce heatmap with protein return filtered dataframe - NSCLC single data", {
  
  CRObject <- getparamhm("nsclc-single")
  output <- do.call(heatmapSC,CRObject)
  
  ggsave("output/nsclc-single_heatmap.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","nsclc-single_heatmap.png")
  
  expect_type(output,"list")
  expected.elements = c("plot","data")
  expect_setequal(names(output), expected.elements)
})    

test_that("Produce heatmap with protein return filtered dataframe - NSCLC multi data", {
  
  CRObject <- getparamhm("nsclc-multi")
  output <- do.call(heatmapSC,CRObject)
  
  ggsave("output/nsclc-multi_heatmap.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","nsclc-single_heatmap.png")
  
  expect_type(output,"list")
  expected.elements = c("plot","data")
  expect_setequal(names(output), expected.elements)
})   

test_that("Produce heatmap with filtered dataframe - BRCA data", {
  
  CRObject <- getparamhm("BRCA")
  output <- do.call(heatmapSC,CRObject)
  
  ggsave("output/nsclc-multi_heatmap.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","nsclc-single_heatmap.png")
  
  expect_type(output,"list")
  expected.elements = c("plot","data")
  expect_setequal(names(output), expected.elements)
})  