test_that("Violin plot returns ggplot object", {
  
  so_demo <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  genes_demo <- read.csv(test_path("fixtures", "Marker_Table_demo.csv"))[,1][1:4]
  
  violin_res <- ViolinPlot(so = so_demo, ident_of_interest = "orig_ident", 
             groups_of_interest = c("1_E13","2_E15","3_Newborn","4_Adult"),
             genes_of_interest = genes_demo)
  
  expected_elements = c("gg","ggplot")
  expect_setequal(class(violin_res), expected_elements)
})
