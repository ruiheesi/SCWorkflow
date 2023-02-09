test_that("Violin plot works for TEC data", {
  
  TEC_data <- select_dataset_SCviolin("TEC")
  
  violin_test <- violinPlot(so = TEC_data$object, ident.of.interest = TEC_data$ident.of.interest, 
             groups.of.interest = TEC_data$groups.of.interest,
             genes.of.interest = TEC_data$genes.of.interest)
  
  expected_elements = c("gg","ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for Chariou data", {
  
  Chariou_data <- select_dataset_SCviolin("Chariou")
  
  violin_test <- violinPlot(so = Chariou_data$object, ident.of.interest = Chariou_data$ident.of.interest, 
                            groups.of.interest = Chariou_data$groups.of.interest,
                            genes.of.interest = Chariou_data$genes.of.interest)
  
  expected_elements = c("gg","ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for NSCLC_Single data", {
  
  NSCLC_Single_data <- select_dataset_SCviolin("NSCLC_Single")
  
  violin_test <- violinPlot(so = NSCLC_Single_data$object, ident.of.interest = NSCLC_Single_data$ident.of.interest, 
                            groups.of.interest = NSCLC_Single_data$groups.of.interest,
                            genes.of.interest = NSCLC_Single_data$genes.of.interest)
  
  expected_elements = c("gg","ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for NSCLC_Multi data", {
  
  NSCLC_Multi_data <- select_dataset_SCviolin("NSCLC_Multi")
  
  violin_test <- violinPlot(so = NSCLC_Multi_data$object, ident.of.interest = NSCLC_Multi_data$ident.of.interest, 
                            groups.of.interest = NSCLC_Multi_data$groups.of.interest,
                            genes.of.interest = NSCLC_Multi_data$genes.of.interest)
  
  expected_elements = c("gg","ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for BRCA data", {
  
  BRCA <- select_dataset_SCviolin("BRCA")
  
  violin_test <- violinPlot(so = BRCA$object, ident.of.interest = BRCA$ident.of.interest, 
                            groups.of.interest = BRCA$groups.of.interest,
                            genes.of.interest = BRCA$genes.of.interest)
  
  expected_elements = c("gg","ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

## Check code detects warnings and errors ##

test_that("Violin plot stops when no query genes are found in the data", {
  
  NSCLC_Single_data <- select_dataset_SCviolin("NSCLC_Single")
  
  expect_error(violinPlot(so = NSCLC_Single_data$object, ident.of.interest = NSCLC_Single_data$ident.of.interest, 
                          groups.of.interest = NSCLC_Single_data$groups.of.interest,
                          genes.of.interest = paste("jibberish", 1:5, sep="_")), "No query genes were found in the dataset.")
  
})

test_that("Violin plot stops when ident of interest is not found in seurat", {
  
  NSCLC_Single_data <- select_dataset_SCviolin("NSCLC_Single")
  
  expect_error(violinPlot(so = NSCLC_Single_data$object, ident.of.interest = "jibberish", 
                          groups.of.interest = NSCLC_Single_data$groups.of.interest,
                          genes.of.interest = NSCLC_Single_data$genes.of.interest), "Unable to find ident of interest in metadata.")
  
})

test_that("Violin plot stops when group of interest is empty", {
  
  NSCLC_Single_data <- select_dataset_SCviolin("NSCLC_Single")
  
  expect_error(violinPlot(so = NSCLC_Single_data$object, ident.of.interest = NSCLC_Single_data$ident.of.interest, 
                          groups.of.interest = paste("jibberish", 1:5, sep="_"),
                          genes.of.interest = NSCLC_Single_data$genes.of.interest), "No groups were found in the selected ident.")
  
})

test_that("Violin plot stops when user attempts to rename ident.of.interest as Gene, Expression, or Scaled", {
  
  NSCLC_Single_data <- select_dataset_SCviolin("NSCLC_Single")
  
  expect_error(violinPlot(so = NSCLC_Single_data$object, ident.of.interest = NSCLC_Single_data$ident.of.interest, 
                          groups.of.interest = NSCLC_Single_data$groups.of.interest,
                          genes.of.interest = NSCLC_Single_data$genes.of.interest,
                          rename.ident = "Gene"), "New ident name cannot be one of Gene, Expression, or scaled.")
  
})

