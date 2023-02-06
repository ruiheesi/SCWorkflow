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
      add.gene.or.protein <- TRUE
      transcripts <- sample(rownames(object),10, replace = FALSE) #Transcripts to plot
      rna.annotations <- transcripts[c(1,2)]
      plot.title <- "Heatmap_Chariou_test"
  } else if (data == "nsclc-single"){
      object <- select_srobject("nsclc-single")
      sample.names <- c("NSCLC_Single")
      metadata <- c("HPCA_main","BP_encode_main")
      set.seed(15)
      transcripts <- sample(rownames(object),10, replace = FALSE) #Transcripts to plot
      add.gene.or.protein <- TRUE
      rna.annotations <- transcripts[c(1,2)]
      #protein.annotations <-  this dataset should have protein.
      plot.title <- "Heatmap_Single_NSCLC"
  } else if (data == "nsclc-multi"){
      object <- select_srobject("nsclc-multi")
      sample.names <- c("Donor_1","Donor_2","Donor_3","Donor_4","Donor_5","Donor_6","Donor_7")
      metadata <- c("HPCA_main","BP_encode_main")
      set.seed(15)
      transcripts <- sample(rownames(object),10, replace = FALSE) #Transcripts to plot
      add.gene.or.protein <- TRUE
      rna.annotations <- transcripts[c(1,2)]
      #protein.annotations <-  this dataset should have protein.
      plot.title <- "Heatmap_Multi_NSCLC"
  } else if (data == "BRCA"){
      object <- select_crobject("BRCA")
      sample.names <- c("CID4471","CID4290A","CID44971","CID4040","CID4513","CID4535")
      metadata <- c("orig.ident","SCT_snn_res.0.2")
      set.seed(15)
      transcripts <- sample(rownames(object@assays$SCT@scale.data),10, replace = FALSE) #Transcripts to plot
      add.gene.or.protein <- TRUE
      rna.annotations <- transcripts[c(1,2)]
      #protein.annotations <-  this dataset should have protein.
      plot.title <- "Heatmap_BRCA"
  }
  
  
  
  return(list("object" = object, 
              "sample.names" = sample.names,
              "metadata"= metadata, 
              "transcripts" = transcripts,
              "plot.title" = plot.title,
              "rna.annotations" = rna.annotations))  
}