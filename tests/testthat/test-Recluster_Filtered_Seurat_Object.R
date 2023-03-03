test_that("Test Recluster Filter Seurat Object using TEC (Mouse) dataset.", {
  
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  output <- reclusterFilteredSeuratObject(object = obj,
                                          old.columns.to.save = c("SCT_snn_res.0.2","SCT_snn_res.0.4","SCT_snn_res.0.6","SCT_snn_res.0.8","SCT_snn_res.1","SCT_snn_res.1.2")
                                          )
  
  #  ggsave("output/TEC_output_umap.pdf",output$plot, width = 10, height = 10)
  #  expect_snapshot_file("output","TEC_output_umap.pdf")
  
  ## Compare
  expect_type(output,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements) 
  }
)




test_that("Test Recluster Filter Seurat Object using TEC (Mouse) dataset with UMAP.", {
  
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  output <- reclusterFilteredSeuratObject(object = obj,
                                          old.columns.to.save = c("SCT_snn_res.0.2","SCT_snn_res.0.4","SCT_snn_res.0.6","SCT_snn_res.0.8","SCT_snn_res.1","SCT_snn_res.1.2"),
                                          reduction.type = "umap"
                                          )
  
  #  ggsave("output/TEC_output_umap.pdf",output$plot, width = 10, height = 10)
  #  expect_snapshot_file("output","TEC_output_umap.pdf")
  
  ## Compare
  expect_type(output,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements) 
  }
)




test_that("Test Recluster Filter Seurat Object using Chariou (Mouse) dataset.", {

  obj <- readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
  output <- reclusterFilteredSeuratObject(object = obj,
                                          old.columns.to.save = c("SCT_snn_res.2.4","SCT_snn_res.2.6","SCT_snn_res.2.8")
                                          )
  
  #  ggsave("output/Chariou_output_umap.pdf",output$plot, width = 10, height = 10)
  #  expect_snapshot_file("output","Chariou_output_umap.pdf")
  
  ## Compare
  expect_type(output,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements) 
  }
)



test_that("Test Recluster Filter Seurat Object using NSCLCSingle (Human) dataset.", {
  
  obj <- readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
  output <- reclusterFilteredSeuratObject(object = obj,
                                          old.columns.to.save = c("SCT_snn_res.0.2","SCT_snn_res.0.4","SCT_snn_res.0.6","SCT_snn_res.0.8","SCT_snn_res.1","SCT_snn_res.1.2")
  )
  
  #  ggsave("output/NSCLCSingle_output_umap.pdf",output$plot, width = 10, height = 10)
  #  expect_snapshot_file("output","NSCLCSingle_output_umap.pdf")
  
  ## Compare
  expect_type(output,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements) 
  }
)



test_that("Test Recluster Filter Seurat Object using NSCLCMulti (Human) dataset.", {
  
  obj <- readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
  output <- reclusterFilteredSeuratObject(object = obj,
                                          old.columns.to.save = c("SCT_snn_res.0.2","SCT_snn_res.0.4","SCT_snn_res.0.6","SCT_snn_res.0.8","SCT_snn_res.1","SCT_snn_res.1.2")
  )
  
  #  ggsave("output/NSCLCMulti_output_umap.pdf",output$plot, width = 10, height = 10)
  #  expect_snapshot_file("output","NSCLCMulti_output_umap.pdf")
  
  ## Compare
  expect_type(output,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements) 
  }
)



test_that("Test Recluster Filter Seurat Object using BRCA (Human) dataset.", {
  
  obj <- readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
  output <- reclusterFilteredSeuratObject(object = obj,
                                          old.columns.to.save = c("SCT_snn_res.0.2","SCT_snn_res.0.4","SCT_snn_res.0.6","SCT_snn_res.0.8","SCT_snn_res.1","SCT_snn_res.1.2")
  )
  
  #  ggsave("output/BRCA_output_umap.pdf",output$plot, width = 10, height = 10)
  #  expect_snapshot_file("output","BRCA_output_umap.pdf")
  
  ## Compare
  expect_type(output,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements) 
  }
)


