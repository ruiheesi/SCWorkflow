
test_that("3D coordinates embedded dataset are returned", {
  
  # load data
  seurat.object <-
    readRDS(test_path("fixtures", "test_tSNE_3D_Coordinates_input.rds"))
 
  # run function
  tSNE_3D_Coordinates <-
    tSNE_3D_Coordinates(Combine_and_Renormalize = seurat.object,
                        Max_sample = 1500)
  ground_truth <- 
    readRDS(test_path("fixtures", "test_tSNE_3D_Coordinates_ground_truth.rds"))
  
  comp_res <- tSNE_3D_Coordinates==ground_truth
  comp_res_nonxyz <- subset(comp_res, select=-c(x,y,z))
  
  # test
  expect_true(sum(comp_res_nonxyz=="FALSE") == 0)
  
})