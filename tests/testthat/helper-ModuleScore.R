
#This helper script will return parameters for each test dependent on the data input:

getparam_moduleScore <- function(data){

  if(data == "TEC"){
    
      object = readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
      sample.names = unique(object$orig.ident)
      sample.to.display = unique(object$orig.ident)
      geneset.dataframe = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      celltypes.to.analyze = colnames(geneset.dataframe)[1:3]
      general.class = celltypes.to.analyze
      levels.df.demo = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))

  } else if (data == "Chariou") {
    
      object = readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
      sample.names = unique(object$orig.ident)
      sample.to.display = unique(object$orig.ident)
      geneset.dataframe = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      celltypes.to.analyze = colnames(geneset.dataframe)[1:3]
      general.class = celltypes.to.analyze
      levels.df.demo = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))
      
  } else if (data == "NSCLC_Single") {
    
    object = readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
    sample.names = unique(object$orig.ident)
    sample.to.display = unique(object$orig.ident)
    geneset.dataframe = data.frame(rand_type1 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type2 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type3 = sample(rownames(object), 5, replace = FALSE))
    celltypes.to.analyze = colnames(geneset.dataframe)
    general.class = celltypes.to.analyze
    levels.df.demo = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))
    
  } else if (data == "NSCLC_Multi") {

    object = readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
    sample.names = unique(object$orig.ident)
    sample.to.display = unique(object$orig.ident)
    geneset.dataframe = data.frame(rand_type1 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type2 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type3 = sample(rownames(object), 5, replace = FALSE))
    celltypes.to.analyze = colnames(geneset.dataframe)
    general.class = celltypes.to.analyze
    levels.df.demo = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))

  } else if (data == "BRCA") {

    object = readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
    sample.names = unique(object$orig.ident)
    sample.to.display = unique(object$orig.ident)
    geneset.dataframe = data.frame(rand_type1 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type2 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type3 = sample(rownames(object), 5, replace = FALSE))
    celltypes.to.analyze = colnames(geneset.dataframe)
    general.class = celltypes.to.analyze
    levels.df.demo = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))
    
  }
  
  return(list("object" = object, "sample.names"= sample.names, 
              "sample.to.display" = sample.to.display, 
              "geneset.dataframe" = geneset.dataframe, 
              "celltypes.to.analyze" = celltypes.to.analyze,
              "general.class" = general.class,
              "levels.df.demo" = levels.df.demo
              ))  
}
