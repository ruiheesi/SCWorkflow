TEC_data <- getparamdp("TEC")

test_that("dotplot run with normal parameters - TEC Data", {
  
  output <- do.call(dotPlotMet,TEC_data)

  ggsave("output/TEC_dotplot.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","TEC_dotplot.png")
  
  expect_type(output,"list")
  expected.elements <- c("plot","pct","exp")
  expect_setequal(names(output), expected.elements)
})  

test_that("dotplot run with warning for fewer categories - TEC Data", {
  
  TEC_sample <- TEC_data
  TEC_sample$cells <- c(0,1,2,3,4,5,6)
  expect_warning(do.call(dotPlotMet,TEC_sample),
                 "^There are")
})

test_that("dotplot run with error for non-matched category - TEC Data", {
  
  TEC_sample <- TEC_data
  TEC_sample$metadata <- "orig.ident"
  expect_error(do.call(dotPlotMet,TEC_sample), 
               "^At least 2 metadata categories you wish to plot should")
})

test_that("dotplot run with warning for duplicate genes", {
  
  TEC_sample <- TEC_data 
  TEC_sample$markers <- c("Map7",TEC_sample$markers)
  expect_warning(do.call(dotPlotMet,TEC_sample), fixed=TRUE, 
                 "The following duplicate genes were removed: Map7")
  
})

test_that("dotplot run with warning for missing genes", {
  
  TEC_sample <- TEC_data 
  TEC_sample$markers <- c("Adgre1", "Ccr2",TEC_sample$markers) 
  expect_warning(do.call(dotPlotMet,TEC_sample), "^There are.")
  
  #Full message:
  #There are 2 gene(s) absent from dataset:'Adgre1', 'Ccr2'. 
  #Possible reasons are that gene is not official gene symbol or 
  #gene is not highly expressed and has been filtered.
})

test_that("dotplot run with error for missing all genes", {
  
  TEC_sample <- TEC_data 
  TEC_sample$markers <- c("APOE","ARG1","CD38","CD3D","CD3E","CD3G")
  expect_error(do.call(dotPlotMet,TEC_sample),
               "No genes listed are found in dataset.")
  
})

test_that("dotplot produced and contingency table returned - Chariou Data", {
  
  Chariou_data <- getparamdp("Chariou")
  output <- do.call(dotPlotMet,Chariou_data)
  
  ggsave("output/Chariou_dotplot.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","Chariou_dotplot.png")
  
  expect_type(output,"list")
  expected.elements <- c("plot","pct","exp")
  expect_setequal(names(output), expected.elements)
  
})  

test_that("dotplot produced and contingency table returned - NSCLC-single Data", {
  
  nsclc_single_data <- getparamdp("nsclc-single")
  output <- do.call(dotPlotMet,nsclc_single_data)
  
  ggsave("output/nsclc-single_dotplot.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","nsclc-single_dotplot.png")
  
  expect_type(output,"list")
  expected.elements <- c("plot","pct","exp")
  expect_setequal(names(output), expected.elements)
  
})  

test_that("dotplot produced and contingency table returned - NSCLC-multi Data", {
  
  nsclc_multi_data <- getparamdp("nsclc-multi")
  expect_warning(output <- do.call(dotPlotMet,nsclc_multi_data),
      "^There are")
  
  ggsave("output/nsclc-multi_dotplot.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","nsclc-multi_dotplot.png")
  
  expect_type(output,"list")
  expected.elements <- c("plot","pct","exp")
  expect_setequal(names(output), expected.elements)
  
})

test_that("dotplot produced and contingency table returned - BRCA Data", {
  
  BRCA_data <- getparamdp("BRCA")
  expect_warning(output <- do.call(dotPlotMet,BRCA_data),
                 "^There are")
  
  ggsave("output/brca_dotplot.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","brca.png")
  
  expect_type(output,"list")
  expected.elements <- c("plot","pct","exp")
  expect_setequal(names(output), expected.elements)
  
})
