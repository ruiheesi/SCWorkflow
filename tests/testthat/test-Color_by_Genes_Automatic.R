test_that("Color by Genes Automatic returns correct class of figures", {
  
  so_demo <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  marker_tab_demo <- read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
  
  Cbg_demo <- color_by_genes(SO = so_demo, 
                             samples_to_include = c("1_E13","2_E15","3_Newborn","4_Adult"), 
                             samples_to_display = c("1_E13","2_E15","3_Newborn","4_Adult"), 
                             marker_list = marker_tab_demo,
                             cells_of_interest = c("CD8_T","CD4_T","Tregs","Macrophages"))
  
  expected_elements <- c("gtable", "gTree", "grob", "gDesc")
  expect_setequal(class(Cbg_demo), expected_elements)
})
