test_that("Load testing dataset", {
  # Seurat_Object <- readRDS('./otherData/PCAnorm.rds')
  
  Seurat_Object <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  
  
  
  Combine_and_Renormalize_out <- Combine_and_Renormalize(Seurat_Object,
                                                     npcs = 15,
                                                     vars_to_regress = c("percent.mt"),
                                                     integratedata = FALSE,
                                                     clust_res_low=0.2,
                                                     clust_res_high = 1.2,
                                                     clust_res_bin = 0.2,
                                                     only_var_genes = FALSE, 
                                                     draw_umap = TRUE,
                                                     draw_tsne = TRUE,
                                                     imageType = "png",
                                                     nfeatures = 2000,
                                                     low_cut = 0.1,
                                                     high_cut = 8,
                                                     low_cut_disp = 1,
                                                     high_cut_disp = 100000,
                                                     selection_method = "vst",
                                                     cell_hashing_data = FALSE,
                                                     project_name = "scRNAProject",
                                                     doMergeData = TRUE,
                                                     seed_for_PCA = 42,
                                                     seed_for_TSNE = 1,
                                                     seed_for_UMAP = 42,
                                                     Do_SCTransform = TRUE,
                                                     Exclude_sample = 0
                                                     
                                 )

  
  

  
  # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_CombNorm.pdf"),
  #        Combine_and_Renormalize_out$plot)
  # plot(Combine_and_Renormalize_out$plot)
  
  # create output
  expected.elements = c("so","plot")
  expect_setequal(names(Combine_and_Renormalize_out), expected.elements)

  # figure slot is a grob
  expect_equal(class(Combine_and_Renormalize_out$plot)[3], 'grob')
  # SO slot contains data
  expect( nrow(Combine_and_Renormalize_out$so@assays$RNA@counts),'> 0' )
  # plot slot contains data
  expect( object.size(Combine_and_Renormalize_out$plot),'> 0' )
  
})


# library(devtools)
# document()
# load_all()
# test_active_file()
# check()
