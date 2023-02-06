# doMergeData doesn't do anything last used in vesion 34 so removed as parameter



for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Multi'
  
  test_that(paste0("Test Combine & Renormalize - Standard (",data," dataset)"), {
    
    
    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_PCA_Norm_SO_downsample.rds")))
  Combine_and_Renormalize_out <- Combine_and_Renormalize(object$so,
                                                     npcs = 15,
                                                     vars.to.regress = c("percent.mt"),
                                                     integratedata = FALSE,
                                                     clust.res.low=0.2,
                                                     clust.res.high = 1.2,
                                                     clust.res.bin = 0.2,
                                                     only.var.genes = FALSE, 
                                                     draw.umap = TRUE,
                                                     draw.tsne = TRUE,
                                                     imageType = "png",
                                                     nfeatures = 2000,
                                                     low.cut = 0.1,
                                                     high.cut = 8,
                                                     low.cut.disp = 1,
                                                     high.cut.disp = 100000,
                                                     selection.method = "vst",
                                                     cell.hashing.data = FALSE,
                                                     project.name = "scRNAProject",
                                                     # doMergeData = TRUE,
                                                     seed.for.PCA = 42,
                                                     seed.for.TSNE = 1,
                                                     seed.for.UMAP = 42,
                                                     SCTransform = TRUE,
                                                     exclude.sample = ""
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

}
 


for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Multi'
  
  test_that(paste0("Test Combine & Renormalize - NoRegress (",data," dataset)"), {
    
    
    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_PCA_Norm_SO_downsample.rds")))
    Combine_and_Renormalize_out <- Combine_and_Renormalize(object$so,
                                                           npcs = 15,
                                                           vars.to.regress = NULL,
                                                           integratedata = FALSE,
                                                           clust.res.low=0.2,
                                                           clust.res.high = 1.2,
                                                           clust.res.bin = 0.2,
                                                           only.var.genes = FALSE, 
                                                           draw.umap = TRUE,
                                                           draw.tsne = TRUE,
                                                           imageType = "png",
                                                           nfeatures = 2000,
                                                           low.cut = 0.1,
                                                           high.cut = 8,
                                                           low.cut_disp = 1,
                                                           high.cut_disp = 100000,
                                                           selection.method = "vst",
                                                           cell.hashing.data = FALSE,
                                                           project.name = "scRNAProject",
                                                           # doMergeData = TRUE,
                                                           seed.for.PCA = 42,
                                                           seed.for.TSNE = 1,
                                                           seed.for.UMAP = 42,
                                                           SCTransform = TRUE,
                                                           exclude.sample = ""
    )

    
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
  
}


for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Multi'
  
  test_that(paste0("Test Combine & Renormalize - only.var.genes (",data," dataset)"), {
    
    
    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_PCA_Norm_SO_downsample.rds")))
    Combine_and_Renormalize_out <- Combine_and_Renormalize(object$so,
                                                           npcs = 15,
                                                           vars.to.regress = NULL,
                                                           integratedata = FALSE,
                                                           clust.res.low=0.2,
                                                           clust.res.high = 1.2,
                                                           clust.res.bin = 0.2,
                                                           only.var.genes = T, 
                                                           draw.umap = TRUE,
                                                           draw.tsne = TRUE,
                                                           imageType = "png",
                                                           nfeatures = 2000,
                                                           low.cut = 0.1,
                                                           high.cut = 8,
                                                           low.cut_disp = 1,
                                                           high.cut_disp = 100000,
                                                           selection.method = "vst",
                                                           cell.hashing.data = FALSE,
                                                           project.name = "scRNAProject",
                                                           # doMergeData = TRUE,
                                                           seed.for.PCA = 42,
                                                           seed.for.TSNE = 1,
                                                           seed.for.UMAP = 42,
                                                           SCTransform = TRUE,
                                                           exclude.sample = ""
    )

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
  
} 
  


for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Multi'
  
  test_that(paste0("Test Combine & Renormalize - exclude.sample (",data," dataset)"), {
    
    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_PCA_Norm_SO_downsample.rds")))
    if (length(names(object$so))>1) {
      sample=object$so[1]%>%names
    }else{ 
      sample=object$so[1]%>%names
    }
   
    
        Combine_and_Renormalize_out <- Combine_and_Renormalize(object$so,
                                                           npcs = 15,
                                                           vars.to.regress = NULL,
                                                           integratedata = FALSE,
                                                           clust.res.low=0.2,
                                                           clust.res.high = 1.2,
                                                           clust.res.bin = 0.2,
                                                           only.var.genes = FALSE, 
                                                           draw.umap = TRUE,
                                                           draw.tsne = TRUE,
                                                           imageType = "png",
                                                           nfeatures = 2000,
                                                           low.cut = 0.1,
                                                           high.cut = 8,
                                                           low.cut_disp = 1,
                                                           high.cut_disp = 100000,
                                                           selection.method = "vst",
                                                           cell.hashing.data = FALSE,
                                                           project.name = "scRNAProject",
                                                           # doMergeData = TRUE,
                                                           seed.for.PCA = 42,
                                                           seed.for.TSNE = 1,
                                                           seed.for.UMAP = 42,
                                                           SCTransform = TRUE,
                                                           exclude.sample = sample
    )
    
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
  
}


for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Multi'
  
  test_that(paste0("Test Combine & Renormalize - selection.method = mean.var.plot (",data," dataset)"), {
    
    
    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_PCA_Norm_SO_downsample.rds")))
    Combine_and_Renormalize_out <- Combine_and_Renormalize(object$so,
                                                           npcs = 15,
                                                           vars.to.regress = NULL,
                                                           integratedata = FALSE,
                                                           clust.res.low=0.2,
                                                           clust.res.high = 1.2,
                                                           clust.res.bin = 0.2,
                                                           only.var.genes = FALSE, 
                                                           draw.umap = TRUE,
                                                           draw.tsne = TRUE,
                                                           imageType = "png",
                                                           nfeatures = 2000,
                                                           low.cut = 0.1,
                                                           high.cut = 8,
                                                           low.cut_disp = 1,
                                                           high.cut_disp = 100000,
                                                           selection.method = "mean.var.plot",
                                                           cell.hashing.data = FALSE,
                                                           project.name = "scRNAProject",
                                                           # doMergeData = TRUE,
                                                           seed.for.PCA = 42,
                                                           seed.for.TSNE = 1,
                                                           seed.for.UMAP = 42,
                                                           SCTransform = TRUE,
                                                           exclude.sample = ""
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
  
}



for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Multi'
  
  test_that(paste0("Test Combine & Renormalize - integratedata (",data," dataset)"), {
    
    
    object <- readRDS(test_path(paste0("fixtures/",data), paste0(data,"_Filtered_PCA_Norm_SO_downsample.rds")))
    Combine_and_Renormalize_out <- Combine_and_Renormalize(object$so,
                                                           npcs = 15,
                                                           vars.to.regress = NULL,
                                                           integratedata = T,
                                                           clust.res.low=0.2,
                                                           clust.res.high = 1.2,
                                                           clust.res.bin = 0.2,
                                                           only.var.genes = FALSE, 
                                                           draw.umap = TRUE,
                                                           draw.tsne = TRUE,
                                                           imageType = "png",
                                                           nfeatures = 2000,
                                                           low.cut = 0.1,
                                                           high.cut = 8,
                                                           low.cut_disp = 1,
                                                           high.cut_disp = 100000,
                                                           selection.method = "vst",
                                                           cell.hashing.data = FALSE,
                                                           project.name = "scRNAProject",
                                                           # doMergeData = T,
                                                           seed.for.PCA = 42,
                                                           seed.for.TSNE = 1,
                                                           seed.for.UMAP = 42,
                                                           SCTransform = TRUE,
                                                           exclude.sample = ""
    )
    
    
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
  
}



# library(devtools)
# document()
# load_all()
# test_active_file()
# check()



# object=object$so
# npcs = 15
# vars.to.regress = NULL
# integratedata = T
# clust.res.low=0.2
# clust.res.high = 1.2
# clust.res.bin = 0.2
# only.var.genes = FALSE
# draw.umap = TRUE
# draw.tsne = TRUE
# imageType = "png"
# nfeatures = 2000
# low.cut = 0.1
# high.cut = 8
# low.cut_disp = 1
# high.cut_disp = 100000
# selection.method = "vst"
# cell.hashing.data = FALSE
# project.name = "scRNAProject"
# # doMergeData = T,
# seed.for.PCA = 42
# seed.for.TSNE = 1
# seed.for.UMAP = 42
# SCTransform = TRUE
# exclude.sample = ""