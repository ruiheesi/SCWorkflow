getParamRaw <- function(data){
  
  if(data == "TEC"){
    
    data.table(Sample_Name=c('5_ABSC_E13_5_TEC','6_ABSC_Newborn_cTEC',
                              '7_ABSC_Adult_cTEC','8_ABSC_E15_cTEC'),
               Gender=c('M','F','N','F'),
               Rename=c("1_Embryo_13_5","3_Newborn","4_Adult","2_Embryo_15"))%>%
      write.table(
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt")),
      sep = '\t')
    
    
    input=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Mouse"
    sample.metadata.table=
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt"))
    sample.name.column='Sample_Name'
    rename.col="Rename" 
    

  
    } else if (data == "Chariou") {
      
      data.table(Sample_Name=c('SCAF1713_1_1','SCAF1714_2_1','SCAF1715_3_1',
                                'SCAF1716_4_1','SCAF1717_5_1'),
                 Rename=c("PBS","ENT","NHSIL12","Combo","CD8dep"))%>%
        write.table(
          test_path(paste0("fixtures/",data,"/",data,"_metadata.txt")),
          sep = '\t')
      
    input=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Mouse"
    sample.metadata.table=
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt"))
    sample.name.column='Sample_Name'
    rename.col="Rename" 
    
    
    
  } else if (data == "NSCLC_Single") {
    
    data.table(Sample_Name=c('PBMC_20k_3p_HT_nextgem_Chromium_X'),
               Rename=c("PBMC_Single"))%>%
      write.table(
        test_path(paste0("fixtures/",data,"/",data,"_metadata.txt")),
        sep = '\t')
    
    input=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Human"
    sample.metadata.table=
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt"))
    sample.name.column='Sample_Name'
    rename.col="Rename" 

    
    
    
  } else if (data == "NSCLC_Multi") {
  
      data.table(Sample_Name=
        c('NSCLC_40k_DTC_3p_HT_nextgem_donor_1_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_2_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_3_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_4_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_5_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_6_count_sample_feature_bc_matrix',
          'NSCLC_40k_DTC_3p_HT_nextgem_donor_7_count_sample_feature_bc_matrix'),
                 Rename=c("Donor_1","Donor_2",
                          "Donor_3","Donor_4",
                          "Donor_5","Donor_6","Donor_7"))%>%
      write.table(
        test_path(paste0("fixtures/",data,"/",data,"_metadata.txt")),
        sep = '\t')      
    input=list.files(  
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Human"
    sample.metadata.table=
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt"))
    sample.name.column='Sample_Name'
    rename.col="Rename" 
    
    
    
    
  } else if (data == "BRCA") {
    input=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Human"
    sample.metadata.table=
      test_path(paste0("fixtures/",data,"/",data,"_metadata.txt"))
    sample.name.column=NULL
    rename.col=NULL
    
    
  }
  
  return(list("input" = input, 
              "organism"=organism, 
              "sample.metadata.table"=sample.metadata.table,
              'sample.name.column'=sample.name.column,
              "rename.col"=rename.col
              ))  
}
