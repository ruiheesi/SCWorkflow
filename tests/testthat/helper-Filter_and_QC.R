getParamFQ <- function(data){
  
  if(data == "TEC"){
    local.file.paths=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Mouse"
    rename = TRUE
    new.sample.names=c("1_Embryo_13_5","3_Newborn","4_Adult","2_Embryo_15")
    mad.mitoch.value=3
    mad.mitoch <- TRUE
    max.mitoch = 10
    protein <- FALSE 
    
    
    
  } else if (data == "Chariou") {
    local.file.paths=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Mouse"
    rename = TRUE
    new.sample.names=c("PBS","ENT","NHSIL12","Combo","CD8dep")
    mad.mitoch.value=3
    mad.mitoch <- FALSE
    max.mitoch = 25
    protein <- FALSE 
    
    
    
  } else if (data == "pbmc-single") {
    local.file.paths=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Human"
    rename = TRUE
    new.sample.names=c("NSCLC_Single")
    mad.mitoch.value=3
    mad.mitoch <- TRUE
    max.mitoch = 25
    protein <- TRUE 
    
    
    
  } else if (data == "nsclc-multi") {
    local.file.paths=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Human"
    rename = TRUE
    new.sample.names=c("Donor_1","Donor_2",
                       "Donor_3","Donor_4",
                       "Donor_5","Donor_6","Donor_7")
    mad.mitoch.value=3
    mad.mitoch <- TRUE
    max.mitoch = 25
    protein <- TRUE 
    
    
    
    
  } else if (data == "BRCA") {
    local.file.paths=list.files(
      test_path(paste0("fixtures/",data,"/h5files")),".h5",full.names = T)
    organism = "Human"
    rename = FALSE
    new.sample.names=c()
    mad.mitoch.value=3
    mad.mitoch <- TRUE
    max.mitoch = 25
    protein <- FALSE 
    
    
  }
  
  return(list("object" = local.file.paths, "organism"=organism, 
              "rename" = rename, "new.sample.names" = new.sample.names, 
              "mad.mitoch.value" = mad.mitoch.value,
              "mad.mitoch" = mad.mitoch,
              "max.mitoch" = max.mitoch,
              "protein" = protein))  
}