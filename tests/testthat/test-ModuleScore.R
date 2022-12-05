test_that("ModuleScore returns metadata with scores and cell calls", {
  
  so_demo <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  marker_tab_demo <- read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
  levels_df_demo = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))
  
  modscore_demo <- ModuleScore(SO = so_demo, 
              sample_names = c("1_E13","2_E15","3_Newborn","4_Adult"),
              sample_to_display <- c("1_E13","2_E15","3_Newborn","4_Adult"),
              geneset_dataframe = marker_tab_demo,
              proteins_presence <- FALSE,
              celltypes_to_analyze <- c("CD8_T","CD4_T","Tregs","Macrophages","M1","M2","Neutrophils","Monocytes","cDCs","pDC","B_cells","NKs"),
              manual_threshold <- c(0),
              general_class <- c("CD8_T","CD4_T","Tregs","Macrophages","M1","M2","Neutrophils","Monocytes","cDCs","pDC","B_cells","NKs"),
              multi_level_class <- FALSE,
              levels_dataframe = levels_df_demo,
              nbins = 10)
  
  expected_elements <- c("Likely_CellType","CD8_T","CD4_T","Tregs","Macrophages","M1","Neutrophils","Monocytes","cDCs","pDC")
  expect_setequal(colnames(modscore_demo[[2]]@meta.data)[30:39], expected_elements)
})
