getModuleScoreParam <- function(data){

  if(data == "tec"){
    
      object = selectCRObject("TEC")
      samples.subset = unique(object$orig.ident)
      sample.to.display = unique(object$orig.ident)
      marker.table = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      celltypes = colnames(marker.table)[1:3]
      general.class = celltypes
      lvl.df = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))

  } else if (data == "chariou") {
    
      object = selectCRObject("Chariou")
      samples.subset = unique(object$orig.ident)
      sample.to.display = unique(object$orig.ident)
      marker.table = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      celltypes = colnames(marker.table)[1:3]
      general.class = celltypes
      lvl.df = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))
      
  } else if (data == "pbmc.single") {
    
    object = selectCRObject("pbmc-single")
    samples.subset = unique(object$orig.ident)
    sample.to.display = unique(object$orig.ident)
    set.seed(114)
    marker.table = data.frame(rand_type1 = sample(rownames(object), 5, 
                                                  replace = FALSE),
                              rand_type2 = sample(rownames(object), 5, 
                                                  replace = FALSE),
                              rand_type3 = sample(rownames(object), 5, 
                                                  replace = FALSE))
    celltypes = colnames(marker.table)
    general.class = celltypes
    lvl.df = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))
    
  } else if (data == "nsclc.multi") {

    object = selectCRObject("nsclc-multi")
    samples.subset = unique(object$orig.ident)
    sample.to.display = unique(object$orig.ident)
    set.seed(214)
    marker.table = data.frame(rand_type1 = sample(rownames(object), 5, 
                                                  replace = FALSE),
                              rand_type2 = sample(rownames(object), 5,
                                                  replace = FALSE),
                              rand_type3 = sample(rownames(object), 5, 
                                                  replace = FALSE))
    celltypes = colnames(marker.table)
    general.class = celltypes
    lvl.df = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))

  } else if (data == "brca") {

    object = selectCRObject("BRCA")
    samples.subset = unique(object$orig.ident)
    sample.to.display = unique(object$orig.ident)
    set.seed(314)
    marker.table = data.frame(rand_type1 = sample(rownames(object), 5, 
                                                  replace = FALSE),
                              rand_type2 = sample(rownames(object), 5, 
                                                  replace = FALSE),
                              rand_type3 = sample(rownames(object), 5, 
                                                  replace = FALSE))
    celltypes = colnames(marker.table)
    general.class = celltypes
    lvl.df = read.csv(test_path("fixtures", "MS_Levels_demo.csv"))
    
  }
  
  return(list("object" = object, 
              "samples.subset"= samples.subset, 
              "sample.to.display" = sample.to.display, 
              "marker.table" = marker.table, 
              "celltypes" = celltypes,
              "general.class" = general.class,
              "lvl.df" = lvl.df
              ))  
}

.drawMSfig <- function(x, width = 10, height = 10){
  path <- tempfile(fileext = ".png")
  ggsave(path, x, width = 10, height = 10)
  print(path)
}
