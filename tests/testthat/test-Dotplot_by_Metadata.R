TEC_data <- getparam("TEC")

test_that("dotplot run with normal parameters - TEC Data", {
  
  results.list <- do.call(DotplotMet,TEC_data)
  expected.elements <- c("plot","pct","exp")
  expect_setequal(names(results.list), expected.elements)
  
})  

test_that("dotplot run with warning for fewer categories - TEC Data", {

  TEC_sample <- TEC_data
  TEC_sample$cells <- c(0,1,2,3,4,5,6)
  expect_warning(do.call(DotplotMet,TEC_sample),
    "^There are")
})

test_that("dotplot run with error for non-matched category - TEC Data", {

  TEC_sample <- TEC_data
  TEC_sample$metadata <- "orig.ident"
  expect_error(do.call(DotplotMet,TEC_sample), 
    "^At least 2 metadata categories you wish to plot should")
})

test_that("dotplot run with warning for duplicate genes", {

  TEC_sample <- TEC_data 
  TEC_sample$markers <- c("Map7",TEC_sample$markers)
  expect_warning(do.call(DotplotMet,TEC_sample), fixed=TRUE, 
          "The following duplicate genes were removed: Map7")

})

test_that("dotplot run with warning for missing genes", {

  TEC_sample <- TEC_data 
  TEC_sample$markers <- c("Adgre1", "Ccr2",TEC_sample$markers) 
  expect_warning(do.call(DotplotMet,TEC_sample),
      "2 genes are absent from dataset:'Adgre1', 'Ccr2'. Possible reasons are that gene is not official gene symbol or gene is not highly expressed and has been filtered.")

})

test_that("dotplot run with error for missing all genes", {

  TEC_sample <- TEC_data 
  TEC_sample$markers <- c("APOE","ARG1","CD38","CD3D","CD3E","CD3G")
  expect_error(do.call(DotplotMet,TEC_sample),
          "No genes listed are found in dataset.")

})

test_that("dotplot produced and contingency table returned - Chariou Data", {
  
  Chariou_data <- getparam("Chariou")
  results.list <- do.call(DotplotMet,Chariou_data)
  
  expected.elements <- c("plot","pct","exp")
  expect_setequal(names(results.list), expected.elements)
  
})  

