
#This helper script will return parameters for each test dependent on the data input:

getparam_ModuleScore <- function(data){

  if(data == "TEC"){
    
      object = readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
      sample_names = unique(object$orig.ident)
      sample_to_display = unique(object$orig.ident)
      geneset_dataframe = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      celltypes_to_analyze = colnames(geneset_dataframe)[1:3]
      general_class = celltypes_to_analyze
      levels_df_demo = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))

  } else if (data == "Chariou") {
    
      object = readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
      sample_names = unique(object$orig.ident)
      sample_to_display = unique(object$orig.ident)
      geneset_dataframe = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      celltypes_to_analyze = colnames(geneset_dataframe)[1:3]
      general_class = celltypes_to_analyze
      levels_df_demo = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))
      
  } else if (data == "NSCLC_Single") {
    
    object = readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
    sample_names = unique(object$orig.ident)
    sample_to_display = unique(object$orig.ident)
    geneset_dataframe = data.frame(rand_type1 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type2 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type3 = sample(rownames(object), 5, replace = FALSE))
    celltypes_to_analyze = colnames(geneset_dataframe)
    general_class = celltypes_to_analyze
    levels_df_demo = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))
    
  } else if (data == "NSCLC_Multi") {

    object = readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
    sample_names = unique(object$orig.ident)
    sample_to_display = unique(object$orig.ident)
    geneset_dataframe = data.frame(rand_type1 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type2 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type3 = sample(rownames(object), 5, replace = FALSE))
    celltypes_to_analyze = colnames(geneset_dataframe)
    general_class = celltypes_to_analyze
    levels_df_demo = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))

  } else if (data == "BRCA") {

    object = readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
    sample_names = unique(object$orig.ident)
    sample_to_display = unique(object$orig.ident)
    geneset_dataframe = data.frame(rand_type1 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type2 = sample(rownames(object), 5, replace = FALSE),
                                   rand_type3 = sample(rownames(object), 5, replace = FALSE))
    celltypes_to_analyze = colnames(geneset_dataframe)
    general_class = celltypes_to_analyze
    levels_df_demo = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))
    
  }
  
  return(list("object" = object, "sample_names"= sample_names, 
              "sample_to_display" = sample_to_display, 
              "geneset_dataframe" = geneset_dataframe, 
              "celltypes_to_analyze" = celltypes_to_analyze,
              "general_class" = general_class,
              "levels_df_demo" = levels_df_demo
              ))  
}
