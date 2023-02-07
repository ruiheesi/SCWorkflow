<<<<<<< HEAD
for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Multi'
  
  test_that(paste0("Test PCA and Normalization - Standard (",data," dataset)"), {
    
    
    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
  PCA_and_Normalization_out <- PCA_and_Normalization(object$so,
                                                     vars.to.regress = c('percent.mt'),
                                                     vars.to.plot = c('percent.mt','nCount_RNA'),
                                                     npcs = 30,
                                                     nfeatures = 2000,
                                                     low.cut = 1,
                                                     high.cut = 8,
                                                     low.cut.disp = 1,
                                                     high.cut.disp = 100000,
                                                     selection.method = "vst",
                                                     jackstraw = FALSE,
                                                     jackstraw.dims = 5,
                                                     methods.PCA = c("Elbow","Marchenko-Pastur"),
                                                     var.threshold = 0.1,
                                                     imagetype = "png"
                                 )
=======
test_that("Load testing dataset", {
  # Seurat_Object <- readRDS('./otherData/PostFilterQC.rds')
  # Seurat_Object <- readRDS('./fixtures/SO_moduleScore.rds')
  
  Seurat_Object <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  
  
  
  PCA_and_Normalization_out <- PCA_and_Normalization(Seurat_Object,
                                                     vars_to_regress = c('percent.mt'),
                                                     vars_to_plot = c('percent.mt','nCount_RNA'),
                                                     npcs = 30,
                                                     nfeatures = 2000,
                                                     low_cut = 1,
                                                     high_cut = 8,
                                                     low_cut_disp = 1,
                                                     high_cut_disp = 100000,
                                                     selection_method = "vst",
                                                     doJackStraw = FALSE,
                                                     JackStraw_dims = 5,
                                                     methods_PCA = c("Elbow","Marchenko-Pastur"),
                                                     var_threshold = 0.1,
                                                     imageType = "png"
                                 )

  
  
  
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_PCAnorm.pdf"),
         # PCA_and_Normalization_out$plot)
  # plot(PCA_and_Normalization_out$plot)
  
<<<<<<< HEAD
  ### Test parameters
=======
  
  
  ### Test parameters
  
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  # create output
  expected.elements = c("so","plot")
  expect_setequal(names(PCA_and_Normalization_out), expected.elements)
  # SO contains object same length as input
<<<<<<< HEAD
  expect_equal(length(PCA_and_Normalization_out$so),length(object$so)) 
=======
  expect_equal(length(PCA_and_Normalization_out$so),length(Seurat_Object)) 
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  # figure slot is a grob
  expect_equal(class(PCA_and_Normalization_out$plot)[3], 'grob')
  # SO slot contains data
  expect( object.size(PCA_and_Normalization_out$so[[1]]@assays$RNA@counts),'> 0' )
  # plot slot contains data
  expect( object.size(PCA_and_Normalization_out$plot),'> 0' )
  
})

<<<<<<< HEAD
}




################################################################################

for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Single'

  test_that(paste0("Test PCA and Normalization - NoRegress (",data," dataset)"), {


    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
    PCA_and_Normalization_out <- PCA_and_Normalization(object$so,
                                                       vars.to.regress = c(),
                                                       vars.to.plot = c('percent.mt','nCount_RNA'),
                                                       npcs = 30,
                                                       nfeatures = 2000,
                                                       low.cut = 1,
                                                       high.cut = 8,
                                                       low.cut.disp = 1,
                                                       high.cut.disp = 100000,
                                                       selection.method = "vst",
                                                       jackstraw = FALSE,
                                                       jackstraw.dims = 5,
                                                       methods.PCA = c("Elbow","Marchenko-Pastur"),
                                                       var.threshold = 0.1,
                                                       imagetype = "png"
    )
    # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_PCAnorm.pdf"),
    # PCA_and_Normalization_out$plot)
    # plot(PCA_and_Normalization_out$plot)

    ### Test parameters
    # create output
    expected.elements = c("so","plot")
    expect_setequal(names(PCA_and_Normalization_out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(PCA_and_Normalization_out$so),length(object))
    # figure slot is a grob
    expect_equal(class(PCA_and_Normalization_out$plot)[3], 'grob')
    # SO slot contains data
    expect( object.size(PCA_and_Normalization_out$so[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(PCA_and_Normalization_out$plot),'> 0' )

  })

}



################################################################################

for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Single'

  test_that(paste0("Test PCA and Normalization - Jackstraw-elbow (",data," dataset)"), {


    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
    PCA_and_Normalization_out <- PCA_and_Normalization(object$so,
                                                       vars.to.regress = c(),
                                                       vars.to.plot = c('percent.mt','nCount_RNA'),
                                                       npcs = 30,
                                                       nfeatures = 2000,
                                                       low.cut = 1,
                                                       high.cut = 8,
                                                       low.cut.disp = 1,
                                                       high.cut.disp = 100000,
                                                       selection.method = "vst",
                                                       jackstraw = T,
                                                       jackstraw.dims = 5,
                                                       methods.PCA = c("Elbow"),
                                                       var.threshold = 0.1,
                                                       imagetype = "png"
    )
    # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_PCAnorm.pdf"),
    # PCA_and_Normalization_out$plot)
    # plot(PCA_and_Normalization_out$plot)

    ### Test parameters
    # create output
    expected.elements = c("so","plot")
    expect_setequal(names(PCA_and_Normalization_out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(PCA_and_Normalization_out$so),length(object))
    # figure slot is a grob
    expect_equal(class(PCA_and_Normalization_out$plot)[3], 'grob')
    # SO slot contains data
    expect( object.size(PCA_and_Normalization_out$so[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(PCA_and_Normalization_out$plot),'> 0' )

  })

}


################################################################################

for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Single'

  test_that(paste0("Test PCA and Normalization - Jackstraw-elbow (",data," dataset)"), {


    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
    PCA_and_Normalization_out <- PCA_and_Normalization(object$so,
                                                       vars.to.regress = c(),
                                                       vars.to.plot = c('percent.mt','nCount_RNA'),
                                                       npcs = 30,
                                                       nfeatures = 2000,
                                                       low.cut = 1,
                                                       high.cut = 8,
                                                       low.cut.disp = 1,
                                                       high.cut.disp = 100000,
                                                       selection.method = "vst",
                                                       jackstraw = T,
                                                       jackstraw.dims = 5,
                                                       methods.PCA = c("Marchenko-Pastur"),
                                                       var.threshold = 0.1,
                                                       imagetype = "png"
    )
    # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_PCAnorm.pdf"),
    # PCA_and_Normalization_out$plot)
    # plot(PCA_and_Normalization_out$plot)

    ### Test parameters
    # create output
    expected.elements = c("so","plot")
    expect_setequal(names(PCA_and_Normalization_out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(PCA_and_Normalization_out$so),length(object))
    # figure slot is a grob
    expect_equal(class(PCA_and_Normalization_out$plot)[3], 'grob')
    # SO slot contains data
    expect( object.size(PCA_and_Normalization_out$so[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(PCA_and_Normalization_out$plot),'> 0' )

  })

}


################################################################################

for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Single'

  test_that(paste0("Test PCA and Normalization - selection- mean.var.plot (",data," dataset)"), {


    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
    PCA_and_Normalization_out <- PCA_and_Normalization(object$so,
                                                       vars.to.regress = c(),
                                                       vars.to.plot = c('percent.mt','nCount_RNA'),
                                                       npcs = 30,
                                                       nfeatures = 2000,
                                                       low.cut = 1,
                                                       high.cut = 8,
                                                       low.cut.disp = 1,
                                                       high.cut.disp = 100000,
                                                       selection.method = "mean.var.plot",
                                                       jackstraw = FALSE,
                                                       jackstraw.dims = 5,
                                                       methods.PCA = c("Elbow","Marchenko-Pastur"),
                                                       var.threshold = 0.1,
                                                       imagetype = "png"
    )
    # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_PCAnorm.pdf"),
    # PCA_and_Normalization_out$plot)
    # plot(PCA_and_Normalization_out$plot)

    ### Test parameters
    # create output
    expected.elements = c("so","plot")
    expect_setequal(names(PCA_and_Normalization_out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(PCA_and_Normalization_out$so),length(object))
    # figure slot is a grob
    expect_equal(class(PCA_and_Normalization_out$plot)[3], 'grob')
    # SO slot contains data
    expect( object.size(PCA_and_Normalization_out$so[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(PCA_and_Normalization_out$plot),'> 0' )

  })

}




################################################################################

for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Single'

  test_that(paste0("Test PCA and Normalization - selection- dispersion (",data," dataset)"), {


    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
    PCA_and_Normalization_out <- PCA_and_Normalization(object$so,
                                                       vars.to.regress = c(),
                                                       vars.to.plot = c('percent.mt','nCount_RNA'),
                                                       npcs = 30,
                                                       nfeatures = 2000,
                                                       low.cut = 1,
                                                       high.cut = 8,
                                                       low.cut.disp = 1,
                                                       high.cut.disp = 100000,
                                                       selection.method = "dispersion",
                                                       jackstraw = FALSE,
                                                       jackstraw.dims = 5,
                                                       methods.PCA = c("Elbow","Marchenko-Pastur"),
                                                       var.threshold = 0.1,
                                                       imagetype = "png"
    )
    # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_PCAnorm.pdf"),
    # PCA_and_Normalization_out$plot)
    # plot(PCA_and_Normalization_out$plot)

    ### Test parameters
    # create output
    expected.elements = c("so","plot")
    expect_setequal(names(PCA_and_Normalization_out), expected.elements)
    # SO contains object same length as input
    expect_equal(length(PCA_and_Normalization_out$so),length(object))
    # figure slot is a grob
    expect_equal(class(PCA_and_Normalization_out$plot)[3], 'grob')
    # SO slot contains data
    expect( object.size(PCA_and_Normalization_out$so[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(PCA_and_Normalization_out$plot),'> 0' )

  })

}


=======
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de

# library(devtools)
# document()
# load_all()
# test_active_file()

<<<<<<< HEAD
# vars.to.regress = c('percent.mt')
# vars.to.plot = c('percent.mt','nCount_RNA')
# npcs = 30
# nfeatures = 2000
# low.cut = 1
# high.cut = 8
# low.cut.disp = 1
# high.cut.disp = 100000
# selection.method = 'vst'
# jackstraw = FALSE
# jackstraw.dims = 5
# methods.PCA = c('Elbow','Marchenko-Pastur')
# var.threshold = 0.1
# imagetype = 'png'
=======
# vars_to_regress = c('percent.mt')
# vars_to_plot = c('percent.mt','nCount_RNA')
# npcs = 30
# nfeatures = 2000
# low_cut = 1
# high_cut = 8
# low_cut_disp = 1
# high_cut_disp = 100000
# selection_method = 'vst'
# doJackStraw = FALSE
# JackStraw_dims = 5
# methods_PCA = c('Elbow','Marchenko-Pastur')
# var_threshold = 0.1
# imageType = 'png'
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de

# saveRDS(PCA_and_Normalization_out$so,'/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/PCAnorm.rds')

# ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_PCAnormgrob.pdf"), grobs)
