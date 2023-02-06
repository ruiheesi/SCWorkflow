## Need to test HTO/split data data
## should I remove Protein data option

for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Single'

test_that(paste0("Test Filter and QC - Standard (",data," dataset)"), {
  if (data=='TEC') {org='Mouse'
  }else{org='Human'}
  
print(paste0("\n Test Filter and QC - Standard (",data," dataset)"))

  localFilePaths=list.files(test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)#%>%as.list

  filter_qc_out <- Filter_and_QC(localFilePaths,
                                 organism = org,
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
                                 )%>%suppressMessages()%>%suppressWarnings()
                            
  
  # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_FilterQC.pdf"),
         # filter_qc_out$plot)
  # plot(filter_qc_out$plot)

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

}

################################################################

for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Single'
  
  test_that(paste0("Test Filter and QC - plot_histogram (",data," dataset)"), {
    if (data=='TEC') {org='Mouse'
    }else{org='Human'}
    
    print(paste0("\n Test Filter and QC - plot_histogram (",data," dataset)"))
    
    localFilePaths=list.files(test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)#%>%as.list
        filter_qc_out <- Filter_and_QC(localFilePaths,
                                   organism = org,
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
                                   File_Filter_Regex = "",
                                   Keep = TRUE,
                                   Split_H5 = FALSE,
                                   Protein = FALSE,
                                   Cell_hash = FALSE,
                                   imageType = "png",
                                   plot_histogram = T
    )%>%suppressMessages()%>%suppressWarnings()
    
    
    # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_FilterQC.pdf"),
    #        filter_qc_out$plot)
    # plot(filter_qc_out$plot)
    
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
  
}

################################################################

for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Single'
  
  test_that(paste0("Test Filter and QC - plot_rename (",data," dataset)"), {
    if (data=='TEC') {org='Mouse'
    }else{org='Human'}
    
    print(paste0("\n Test Filter and QC - rename (",data," dataset)"))
    
    localFilePaths=list.files(test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)#%>%as.list
    filter_qc_out <- Filter_and_QC(localFilePaths,
                                   organism = org,
                                   rename = T,
                                   New_Sample_Names = paste0('Sample_',c(1:length(localFilePaths))),
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
                                   plot_histogram = F
    )%>%suppressMessages()%>%suppressWarnings()
    
    
    # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_FilterQC.pdf"),
    # filter_qc_out$plot)
    # plot(filter_qc_out$plot)
    
    # create output
    expected.elements = c("so","plot")
    # expect_setequal(names(filter_qc_out), expected.elements)
    expect_setequal(names(filter_qc_out$so), paste0('Sample_',c(1:length(localFilePaths))))
    # SO contains object same length as input
    expect_equal(length(filter_qc_out$so),length(localFilePaths))
    # figure slot is a grob
    expect_equal(class(filter_qc_out$plot)[3], 'grob')
    # SO slot contains data
    expect( object.size(filter_qc_out$so[[1]]@assays$RNA@counts),'> 0' )
    # plot slot contains data
    expect( object.size(filter_qc_out$plot),'> 0' )
    
  })
  
}

################################################################

for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='NSCLC_Single'
  
  test_that(paste0("Test Filter and QC - Filter_VDJ_Genes (",data," dataset)"), {
    if (data=='TEC') {org='Mouse'
    }else{org='Human'}
    
    print(paste0("\n Test Filter and QC - Filter_VDJ_Genes (",data," dataset)"))
    
    localFilePaths=list.files(test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)#%>%as.list
    filter_qc_out <- Filter_and_QC(localFilePaths,
                                   organism = org,
                                   rename = F,
                                   New_Sample_Names = c(),
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
                                   Filter_VDJ_Genes = T,
                                   Keep = TRUE,
                                   File_Filter_Regex = "",
                                   Split_H5 = FALSE,
                                   Protein = FALSE,
                                   Cell_hash = FALSE,
                                   imageType = "png",
                                   plot_histogram = F
    )%>%suppressMessages()%>%suppressWarnings()
    
    
    # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_FilterQC.pdf"),
    # filter_qc_out$plot)
    # plot(filter_qc_out$plot)
    
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
  
}

################################################################

for (data in c('TEC','Chariou','NSCLC_Single','NSCLC_Multi')) {
  # data='TEC'
  
  test_that(paste0("Test Filter and QC - Regex (",data," dataset)"), {
    if (data=='TEC') {org='Mouse'
    }else{org='Human'}
    
    print(paste0("\n Test Filter and QC - Regex (",data," dataset)"))
    
    localFilePaths=list.files(test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)#%>%as.list
    filter_qc_out <- Filter_and_QC(localFilePaths,
                                   organism = org,
                                   rename = T,
                                   New_Sample_Names = paste0('Sample_',c(1:length(localFilePaths))),
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
                                   Filter_VDJ_Genes = F,
                                   Keep = TRUE,
                                   File_Filter_Regex = "Sample_[1-2]",
                                   Split_H5 = FALSE,
                                   Protein = FALSE,
                                   Cell_hash = FALSE,
                                   imageType = "png",
                                   plot_histogram = F
    )%>%suppressMessages()%>%suppressWarnings()
    
    

    
    # ggsave(file=paste0("/rstudio-files/ccbr-data/users/phil/SCWorkflow/tests/testthat/otherData/test_FilterQC.pdf"),
    # filter_qc_out$plot)
    # plot(filter_qc_out$plot)
    
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
  
}


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
