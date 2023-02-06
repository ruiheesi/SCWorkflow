## Check Color by Genes Automatic works for different datasets ##

test_that("Color by Genes Automatic works for TEC data", {

  TEC_data <- getparam_Cbg_auto("TEC")

  Cbg_demo <- color_by_genes(SO = TEC_data$object,
                             samples_to_include = TEC_data$samples_to_include,
                             samples_to_display = TEC_data$samples_to_display,
                             marker_list = TEC_data$marker_list,
                             cells_of_interest = TEC_data$cells_of_interest)

  expected_elements <- c("gtable", "gTree", "grob", "gDesc")
  expect_setequal(class(Cbg_demo), expected_elements)
})

test_that("Color by Genes Automatic works for Chariou data", {
  
  Chariou_data <- getparam_Cbg_auto("Chariou")
  
  Cbg_demo <- color_by_genes(SO = Chariou_data$object,
                             samples_to_include = Chariou_data$samples_to_include,
                             samples_to_display = Chariou_data$samples_to_display,
                             marker_list = Chariou_data$marker_list,
                             cells_of_interest = Chariou_data$cells_of_interest)
  
  expected_elements <- c("gtable", "gTree", "grob", "gDesc")
  expect_setequal(class(Cbg_demo), expected_elements)
  
})

test_that("Color by Genes Automatic works for NSCLC_Single data", {
  
  NSCLC_Single_data <- getparam_Cbg_auto("NSCLC_Single")
  
  Cbg_demo <- color_by_genes(SO = NSCLC_Single_data$object,
                             samples_to_include = NSCLC_Single_data$samples_to_include,
                             samples_to_display = NSCLC_Single_data$samples_to_display,
                             marker_list = NSCLC_Single_data$marker_list,
                             cells_of_interest = NSCLC_Single_data$cells_of_interest)
  
  expected_elements <- c("gtable", "gTree", "grob", "gDesc")
  expect_setequal(class(Cbg_demo), expected_elements)
  
})

test_that("Color by Genes Automatic works for NSCLC_Multi data", {
  
  NSCLC_Multi_data <- getparam_Cbg_auto("NSCLC_Multi")
  
  Cbg_demo <- color_by_genes(SO = NSCLC_Multi_data$object,
                             samples_to_include = NSCLC_Multi_data$samples_to_include,
                             samples_to_display = NSCLC_Multi_data$samples_to_display,
                             marker_list = NSCLC_Multi_data$marker_list,
                             cells_of_interest = NSCLC_Multi_data$cells_of_interest)
  
  expected_elements <- c("gtable", "gTree", "grob", "gDesc")
  expect_setequal(class(Cbg_demo), expected_elements)
  
})

test_that("Color by Genes Automatic works for BRCA data", {
  
  BRCA_data <- getparam_Cbg_auto("BRCA")
  
  Cbg_demo <- color_by_genes(SO = BRCA_data$object,
                             samples_to_include = BRCA_data$samples_to_include,
                             samples_to_display = BRCA_data$samples_to_display,
                             marker_list = BRCA_data$marker_list,
                             cells_of_interest = BRCA_data$cells_of_interest)
  
  expected_elements <- c("gtable", "gTree", "grob", "gDesc")
  expect_setequal(class(Cbg_demo), expected_elements)
  
})


## Error Checking ##
test_that("Color by Genes Automatic stops when user inputs an assay not found in seurat ", {
  
  TEC_data <- getparam_Cbg_auto("TEC")
  
  expect_error(color_by_genes(SO = TEC_data$object,
                             samples_to_include = TEC_data$samples_to_include,
                             samples_to_display = TEC_data$samples_to_display,
                             marker_list = TEC_data$marker_list,
                             cells_of_interest = TEC_data$cells_of_interest,
                             assay = "wrong_assay"), "assay type not found in seurat")

})

test_that("Color by Genes Automatic stops when user inputs an reduction type not found in seurat ", {
  
  TEC_data <- getparam_Cbg_auto("TEC")
  
  expect_error(color_by_genes(SO = TEC_data$object,
                              samples_to_include = TEC_data$samples_to_include,
                              samples_to_display = TEC_data$samples_to_display,
                              marker_list = TEC_data$marker_list,
                              cells_of_interest = TEC_data$cells_of_interest,
                              reduction = "wrong_reduction"), "reduction type not found in seurat")
  
})
