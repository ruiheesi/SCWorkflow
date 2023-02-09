test_that("moduleScore returns metadata with scores and cell calls for TEC", {
  
  TEC <- getparam_moduleScore("TEC")

  suppressWarnings(modscore_demo <- moduleScore(SO = TEC$object, 
              sample.names = TEC$sample.names,
              sample.to.display = TEC$sample.to.display,
              geneset.dataframe = TEC$geneset.dataframe,
              celltypes.to.analyze = TEC$celltypes.to.analyze,
              general.class = TEC$general.class,
              levels.dataframe = TEC$levels_df_demo,
              nbins = 10))
  
  expected_elements <- c("Likely_CellType",TEC$celltypes.to.analyze)
  expect(all(expected_elements %in% colnames(modscore_demo$Seurat_Object@meta.data)), 
                          failure_message = "modscore results not found")
  
})

test_that("moduleScore returns metadata with scores and cell calls for Chariou", {
  
  Chariou <- getparam_moduleScore("Chariou")
  
  suppressWarnings(modscore_demo <- moduleScore(SO = Chariou$object, 
                                                sample.names = Chariou$sample.names,
                                                sample.to.display = Chariou$sample.to.display,
                                                geneset.dataframe = Chariou$geneset.dataframe,
                                                celltypes.to.analyze = Chariou$celltypes.to.analyze,
                                                general.class = Chariou$general.class,
                                                levels.dataframe = Chariou$levels_df_demo,
                                                nbins = 10))
  
  expected_elements <- c("Likely_CellType",Chariou$celltypes.to.analyze)
  expect(all(expected_elements %in% colnames(modscore_demo$Seurat_Object@meta.data)), 
         failure_message = "modscore results not found")
  
})

test_that("moduleScore returns metadata with scores and cell calls for NSCLC_Single", {
  
  NSCLC_Single <- getparam_moduleScore("NSCLC_Single")
  
  suppressWarnings(modscore_demo <- moduleScore(SO = NSCLC_Single$object, 
                                                sample.names = NSCLC_Single$sample.names,
                                                sample.to.display = NSCLC_Single$sample.to.display,
                                                geneset.dataframe = NSCLC_Single$geneset.dataframe,
                                                celltypes.to.analyze = NSCLC_Single$celltypes.to.analyze,
                                                general.class = NSCLC_Single$general.class,
                                                levels.dataframe = NSCLC_Single$levels_df_demo,
                                                nbins = 10))
  
  expected_elements <- c("Likely_CellType",NSCLC_Single$celltypes.to.analyze)
  expect(all(expected_elements %in% colnames(modscore_demo$Seurat_Object@meta.data)), 
         failure_message = "modscore results not found")
  
})

test_that("moduleScore returns metadata with scores and cell calls for NSCLC_Multi", {
  
  NSCLC_Multi <- getparam_moduleScore("NSCLC_Multi")
  
  suppressWarnings(modscore_demo <- moduleScore(SO = NSCLC_Multi$object, 
                                                sample.names = NSCLC_Multi$sample.names,
                                                sample.to.display = NSCLC_Multi$sample.to.display,
                                                geneset.dataframe = NSCLC_Multi$geneset.dataframe,
                                                celltypes.to.analyze = NSCLC_Multi$celltypes.to.analyze,
                                                general.class = NSCLC_Multi$general.class,
                                                levels.dataframe = NSCLC_Multi$levels_df_demo,
                                                nbins = 10))
  
  expected_elements <- c("Likely_CellType",NSCLC_Multi$celltypes.to.analyze)
  expect(all(expected_elements %in% colnames(modscore_demo$Seurat_Object@meta.data)), 
         failure_message = "modscore results not found")
  
})

test_that("moduleScore returns metadata with scores and cell calls for BRCA", {
  
  BRCA <- getparam_moduleScore("BRCA")
  
  suppressWarnings(modscore_demo <- moduleScore(SO = BRCA$object, 
                                                sample.names = BRCA$sample.names,
                                                sample.to.display = BRCA$sample.to.display,
                                                geneset.dataframe = BRCA$geneset.dataframe,
                                                celltypes.to.analyze = BRCA$celltypes.to.analyze,
                                                general.class = BRCA$general.class,
                                                levels.dataframe = BRCA$levels_df_demo,
                                                nbins = 10))
  
  expected_elements <- c("Likely_CellType",BRCA$celltypes.to.analyze)
  expect(all(expected_elements %in% colnames(modscore_demo$Seurat_Object@meta.data)), 
         failure_message = "modscore results not found")
  
})

## Error Testings ##

test_that("moduleScore detects when no genes are found in the data", {
  
  Chariou <- getparam_moduleScore("Chariou")
  
  expect_error(suppressWarnings(modscore_demo <- moduleScore(SO = Chariou$object, 
                                                sample.names = Chariou$sample.names,
                                                sample.to.display = Chariou$sample.to.display,
                                                geneset.dataframe = apply(Chariou$geneset.dataframe,2, function(x) toupper(x)),
                                                celltypes.to.analyze = Chariou$celltypes.to.analyze,
                                                general.class = Chariou$general.class,
                                                levels.dataframe = Chariou$levels_df_demo,
                                                nbins = 10), "No genes from list was found in data"))
  
})

test_that("moduleScore detects when threshold number does not match number of cells to analyze", {
  
  Chariou <- getparam_moduleScore("Chariou")
  
  expect_error(suppressWarnings(modscore_demo <- moduleScore(SO = Chariou$object, 
                                                             sample.names = Chariou$sample.names,
                                                             sample.to.display = Chariou$sample.to.display,
                                                             geneset.dataframe = Chariou$geneset.dataframe,
                                                             celltypes.to.analyze = Chariou$celltypes.to.analyze,
                                                             general.class = Chariou$general.class,
                                                             levels.dataframe = Chariou$levels_df_demo,
                                                             manual.threshold = rep(0.1,5),
                                                             nbins = 10), "Manual threshold length does not match number of celltypes to analyze - please check manual thresholds"))
  
})
