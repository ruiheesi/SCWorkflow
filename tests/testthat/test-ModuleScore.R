test_that("ModuleScore returns metadata with scores and cell calls for TEC", {
  
  TEC <- getparam_ModuleScore("TEC")

  suppressWarnings(modscore_demo <- ModuleScore(SO = TEC$object, 
              sample_names = TEC$sample_names,
              sample_to_display = TEC$sample_to_display,
              geneset_dataframe = TEC$geneset_dataframe,
              celltypes_to_analyze = TEC$celltypes_to_analyze,
              general_class = TEC$general_class,
              levels_dataframe = TEC$levels_df_demo,
              nbins = 10))
  
  expected_elements <- c("Likely_CellType",TEC$celltypes_to_analyze)
  expect(all(expected_elements %in% colnames(modscore_demo$Seurat_Object@meta.data)), 
                          failure_message = "modscore results not found")
  
})

test_that("ModuleScore returns metadata with scores and cell calls for Chariou", {
  
  Chariou <- getparam_ModuleScore("Chariou")
  
  suppressWarnings(modscore_demo <- ModuleScore(SO = Chariou$object, 
                                                sample_names = Chariou$sample_names,
                                                sample_to_display = Chariou$sample_to_display,
                                                geneset_dataframe = Chariou$geneset_dataframe,
                                                celltypes_to_analyze = Chariou$celltypes_to_analyze,
                                                general_class = Chariou$general_class,
                                                levels_dataframe = Chariou$levels_df_demo,
                                                nbins = 10))
  
  expected_elements <- c("Likely_CellType",Chariou$celltypes_to_analyze)
  expect(all(expected_elements %in% colnames(modscore_demo$Seurat_Object@meta.data)), 
         failure_message = "modscore results not found")
  
})

test_that("ModuleScore returns metadata with scores and cell calls for NSCLC_Single", {
  
  NSCLC_Single <- getparam_ModuleScore("NSCLC_Single")
  
  suppressWarnings(modscore_demo <- ModuleScore(SO = NSCLC_Single$object, 
                                                sample_names = NSCLC_Single$sample_names,
                                                sample_to_display = NSCLC_Single$sample_to_display,
                                                geneset_dataframe = NSCLC_Single$geneset_dataframe,
                                                celltypes_to_analyze = NSCLC_Single$celltypes_to_analyze,
                                                general_class = NSCLC_Single$general_class,
                                                levels_dataframe = NSCLC_Single$levels_df_demo,
                                                nbins = 10))
  
  expected_elements <- c("Likely_CellType",NSCLC_Single$celltypes_to_analyze)
  expect(all(expected_elements %in% colnames(modscore_demo$Seurat_Object@meta.data)), 
         failure_message = "modscore results not found")
  
})

test_that("ModuleScore returns metadata with scores and cell calls for NSCLC_Multi", {
  
  NSCLC_Multi <- getparam_ModuleScore("NSCLC_Multi")
  
  suppressWarnings(modscore_demo <- ModuleScore(SO = NSCLC_Multi$object, 
                                                sample_names = NSCLC_Multi$sample_names,
                                                sample_to_display = NSCLC_Multi$sample_to_display,
                                                geneset_dataframe = NSCLC_Multi$geneset_dataframe,
                                                celltypes_to_analyze = NSCLC_Multi$celltypes_to_analyze,
                                                general_class = NSCLC_Multi$general_class,
                                                levels_dataframe = NSCLC_Multi$levels_df_demo,
                                                nbins = 10))
  
  expected_elements <- c("Likely_CellType",NSCLC_Multi$celltypes_to_analyze)
  expect(all(expected_elements %in% colnames(modscore_demo$Seurat_Object@meta.data)), 
         failure_message = "modscore results not found")
  
})

test_that("ModuleScore returns metadata with scores and cell calls for BRCA", {
  
  BRCA <- getparam_ModuleScore("BRCA")
  
  suppressWarnings(modscore_demo <- ModuleScore(SO = BRCA$object, 
                                                sample_names = BRCA$sample_names,
                                                sample_to_display = BRCA$sample_to_display,
                                                geneset_dataframe = BRCA$geneset_dataframe,
                                                celltypes_to_analyze = BRCA$celltypes_to_analyze,
                                                general_class = BRCA$general_class,
                                                levels_dataframe = BRCA$levels_df_demo,
                                                nbins = 10))
  
  expected_elements <- c("Likely_CellType",BRCA$celltypes_to_analyze)
  expect(all(expected_elements %in% colnames(modscore_demo$Seurat_Object@meta.data)), 
         failure_message = "modscore results not found")
  
})

## Error Testings ##

test_that("ModuleScore detects when no genes are found in the data", {
  
  Chariou <- getparam_ModuleScore("Chariou")
  
  expect_error(suppressWarnings(modscore_demo <- ModuleScore(SO = Chariou$object, 
                                                sample_names = Chariou$sample_names,
                                                sample_to_display = Chariou$sample_to_display,
                                                geneset_dataframe = apply(Chariou$geneset_dataframe,2, function(x) toupper(x)),
                                                celltypes_to_analyze = Chariou$celltypes_to_analyze,
                                                general_class = Chariou$general_class,
                                                levels_dataframe = Chariou$levels_df_demo,
                                                nbins = 10), "No genes from list was found in data"))
  
})

test_that("ModuleScore detects when threshold number does not match number of cells to analyze", {
  
  Chariou <- getparam_ModuleScore("Chariou")
  
  expect_error(suppressWarnings(modscore_demo <- ModuleScore(SO = Chariou$object, 
                                                             sample_names = Chariou$sample_names,
                                                             sample_to_display = Chariou$sample_to_display,
                                                             geneset_dataframe = Chariou$geneset_dataframe,
                                                             celltypes_to_analyze = Chariou$celltypes_to_analyze,
                                                             general_class = Chariou$general_class,
                                                             levels_dataframe = Chariou$levels_df_demo,
                                                             manual_threshold = rep(0.1,5),
                                                             nbins = 10), "Manual threshold length does not match number of celltypes to analyze - please check manual thresholds"))
  
})
