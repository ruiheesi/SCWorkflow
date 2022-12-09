test_that("Produce heatmap and return filtered dataframe", {
  
      seurat_object <- readRDS("/rstudio-files/ccbr-data/users/maggie/SCWorkflow/tests/testthat/fixtures/SO_moduleScore.rds")
      sample_names <- c("1_E13","2_E15","3_Newborn","4_Adult")
      metadata_to_plot <- c("orig_ident","Likely_CellType")
      transcripts_to_plot = c("Epcam","Aire","Fezf2","Pigr","Ly6d","Spink5","Ivl","Krt10","Gapdh","Cd8a","Foxp3","Cd4")
      proteins_to_plot = c("")
      plot_title <- "Heatmap_IO_test"
      trim_outliers_percentage <- 0.01
      
      heatplot <- Heatmap(object = seurat_object,
                          sample.names = sample_names,
                          metadata = metadata_to_plot,
                          transcripts = transcripts_to_plot,
                          proteins = NULL,
                          trim.outliers.percentage = trim_outliers_percentage,
                          plot.title = "Heatmap")
      print(heatplot$data[1:5,1:5])
      expected.elements = c("plot","data")
      expect_setequal(names(heatplot), expected.elements)
})      
      