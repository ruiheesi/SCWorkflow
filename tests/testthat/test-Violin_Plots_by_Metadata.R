test_that("Violin plot works for TEC data", {
  
  TEC_data <- select_dataset_SCviolin("TEC")
  
  violin_test <- ViolinPlot(so = TEC_data$object, ident_of_interest = TEC_data$ident_of_interest, 
             groups_of_interest = TEC_data$groups_of_interest,
             genes_of_interest = TEC_data$genes_of_interest)
  
  expected_elements = c("gg","ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for Chariou data", {
  
  Chariou_data <- select_dataset_SCviolin("Chariou")
  
  violin_test <- ViolinPlot(so = Chariou_data$object, ident_of_interest = Chariou_data$ident_of_interest, 
                            groups_of_interest = Chariou_data$groups_of_interest,
                            genes_of_interest = Chariou_data$genes_of_interest)
  
  expected_elements = c("gg","ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for NSCLC_Single data", {
  
  NSCLC_Single_data <- select_dataset_SCviolin("NSCLC_Single")
  
  violin_test <- ViolinPlot(so = NSCLC_Single_data$object, ident_of_interest = NSCLC_Single_data$ident_of_interest, 
                            groups_of_interest = NSCLC_Single_data$groups_of_interest,
                            genes_of_interest = NSCLC_Single_data$genes_of_interest)
  
  expected_elements = c("gg","ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for NSCLC_Multi data", {
  
  NSCLC_Multi_data <- select_dataset_SCviolin("NSCLC_Multi")
  
  violin_test <- ViolinPlot(so = NSCLC_Multi_data$object, ident_of_interest = NSCLC_Multi_data$ident_of_interest, 
                            groups_of_interest = NSCLC_Multi_data$groups_of_interest,
                            genes_of_interest = NSCLC_Multi_data$genes_of_interest)
  
  expected_elements = c("gg","ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for BRCA data", {
  
  BRCA <- select_dataset_SCviolin("BRCA")
  
  violin_test <- ViolinPlot(so = BRCA$object, ident_of_interest = BRCA$ident_of_interest, 
                            groups_of_interest = BRCA$groups_of_interest,
                            genes_of_interest = BRCA$genes_of_interest)
  
  expected_elements = c("gg","ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

## Check code detects warnings and errors ##

test_that("Violin plot stops when no query genes are found in the data", {
  
  NSCLC_Single_data <- select_dataset_SCviolin("NSCLC_Single")
  
  expect_error(ViolinPlot(so = NSCLC_Single_data$object, ident_of_interest = NSCLC_Single_data$ident_of_interest, 
                          groups_of_interest = NSCLC_Single_data$groups_of_interest,
                          genes_of_interest = paste("jibberish", 1:5, sep="_")), "No query genes were found in the dataset.")
  
})

test_that("Violin plot stops when ident of interest is not found in seurat", {
  
  NSCLC_Single_data <- select_dataset_SCviolin("NSCLC_Single")
  
  expect_error(ViolinPlot(so = NSCLC_Single_data$object, ident_of_interest = "jibberish", 
                          groups_of_interest = NSCLC_Single_data$groups_of_interest,
                          genes_of_interest = NSCLC_Single_data$genes_of_interest), "Unable to find ident of interest in metadata.")
  
})

test_that("Violin plot stops when group of interest is empty", {
  
  NSCLC_Single_data <- select_dataset_SCviolin("NSCLC_Single")
  
  expect_error(ViolinPlot(so = NSCLC_Single_data$object, ident_of_interest = NSCLC_Single_data$ident_of_interest, 
                          groups_of_interest = paste("jibberish", 1:5, sep="_"),
                          genes_of_interest = NSCLC_Single_data$genes_of_interest), "No groups were found in the selected ident.")
  
})

test_that("Violin plot stops when user attempts to rename ident_of_interest as Gene, Expression, or Scaled", {
  
  NSCLC_Single_data <- select_dataset_SCviolin("NSCLC_Single")
  
  expect_error(ViolinPlot(so = NSCLC_Single_data$object, ident_of_interest = NSCLC_Single_data$ident_of_interest, 
                          groups_of_interest = NSCLC_Single_data$groups_of_interest,
                          genes_of_interest = NSCLC_Single_data$genes_of_interest,
                          rename_ident = "Gene"), "New ident name cannot be one of Gene, Expression, or scaled.")
  
})

