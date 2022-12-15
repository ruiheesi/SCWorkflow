test_that("Load testing dataset", {
  datadir <- readRDS(test_path("fixtures","filter_qc_test_in.rds"))
  
  filter_qc_out <- Filter_and_QC(datadir,
                                 Organism = "Mouse")
                            
  expected.elements = c("filter_qc_test_h5")
  expect_setequal(names(filter_qc_out), expected.elements)

})
