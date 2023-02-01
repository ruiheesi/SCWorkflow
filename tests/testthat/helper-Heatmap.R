#This helper script will return parameters for each test dependent on the data input:
getparamhm <- function(data){

  if(data == "TEC"){
      object <- select_crobject("TEC") 
      sample.names <- c("1_Embryo_13_5","2_Embryo_15","3_Newborn","4_Adult")
      metadata <- c("orig.ident", "Phase")
      set.seed(15)
      transcripts <- sample(rownames(object),10, replace = FALSE) #Transcripts to plot
      rna.annotations <- transcripts[c(1,2)]
      plot.title <- "Heatmap_TEC_test"
  } else if (data == "Chariou") {
      object <- select_crobject("Chariou")
      sample.names <- c("PBS","CD8dep","ENT","NHSIL12","Combo")
      metadata <- "orig.ident"
      set.seed(15)
      transcripts <- sample(rownames(object),10, replace = FALSE) #Transcripts to plot
      rna.annotations <- transcripts[c(1,2)]
      plot.title <- "Heatmap_Chariou_test"
  }
  return(list("object" = object, 
              "sample.names" = sample.names,
              "metadata"= metadata, 
              "transcripts" = transcripts,
              "plot.title" = plot.title,
              "rna.annotations" = rna.annotations))  
}