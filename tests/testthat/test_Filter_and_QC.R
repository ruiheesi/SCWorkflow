test_that("Load testing dataset", {

  datadir <- readRDS(test_path("fixtures","filter_qc_test_in.rds"))

  localFilePaths=(test_path("fixtures", "filter_qc_test_h5.h5"))
  
  
  filter_qc_out <- Filter_and_QC(localFilePaths,
                                 organism = "Mouse",
                                 rename = F,
                                 New_Sample_Names = c("Sample_1", "Sample_2"),
                                 mincells = 10,
                                 mingenes = 500,
                                 complexity = 0.6,
                                 MAD_genes_value = 3,
                                 MAD_mitoch_value = 3,
                                 minUMI = 500,
                                 MAD_gene = TRUE,
                                 MAD_mitoch = TRUE,
                                 maxgenes = 2500,
                                 maxmitoch = 10,
                                 Filter_VDJ_Genes = FALSE,
                                 Keep = TRUE,
                                 File_Filter_Regex = "",
                                 Split_H5 = FALSE,
                                 Protein = FALSE,
                                 Cell_hash = FALSE,
                                 imageType = "png",
                                 plot_histogram = FALSE
                                 )

                            
  # create output
  expected.elements = c("so","plot")
  expect_setequal(names(filter_qc_out), expected.elements)
  # SO contains object same length as input
  expect_equal(length(filter_qc_out$so),length(localFilePaths))
  # figure slot is a grob
  expect_equal(class(filter_qc_out$plot)[3], 'grob')
  # SO slot contains data
  expect( object.size(filter_qc_out$so[[1]]@assays$RNA@counts),'> 0' )
  # plot slot contains data
  expect( object.size(filter_qc_out$plot),'> 0' )

})

# library(devtools)
# document()
# load_all()
# test_active_file()

# organism = "Mouse"
# rename = F
# New_Sample_Names = c("Sample_1", "Sample_2")
# mincells = 10
# mingenes = 500
# complexity = 0.6
# MAD_genes_value = 3
# MAD_mitoch_value = 3
# minUMI = 500
# MAD_gene = TRUE
# MAD_mitoch = TRUE
# maxgenes = 2500
# maxmitoch = 10
# Filter_VDJ_Genes = FALSE
# Keep = TRUE
# File_Filter_Regex = ""
# Split_H5 = FALSE
# Protein = FALSE
# Cell_hash = FALSE
# imageType = "png"
# plot_histogram = FALSE

# saveRDS(filter_qc_out$so,'/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/FilterQC.rds')
