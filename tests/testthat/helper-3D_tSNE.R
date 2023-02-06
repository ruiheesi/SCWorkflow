#This helper script will return parameters for each test dependent on the data input:
getparam3d <- function(data){

  if(data == "TEC"){
      object <- select_srobject("TEC")
      color.variable = "orig.ident"
      label.variable = "immgen_main"
      filename = "output/TEC_tsneplot.html"
      save.plot = TRUE
  } else if (data == "Chariou"){
    object <- select_srobject("Chariou")
    color.variable = "mouseRNAseq"
    label.variable = "Barcode"
    filename = "output/Chariou_tsneplot.html"
    save.plot = TRUE
  } else if (data == "nsclc-single"){
    object <- select_srobject("nsclc-single")
    color.variable = "SCT_snn_res.0.2"
    label.variable = "BP_encode_main"
    filename = "output/nsclc-single_tsneplot.html"
    save.plot = TRUE
  } else if (data == "nsclc-multi"){
    object <- select_srobject("nsclc-multi")
    color.variable = "orig.ident"
    label.variable = "BP_encode_main"
    filename = "output/nsclc-multi_tsneplot.html"
    save.plot = TRUE
  } else if (data == "BRCA"){
    object <- select_crobject("BRCA")
    color.variable = "orig.ident"
    label.variable = "SCT_snn_res.0.4"
    filename = "output/BRCA_tsneplot.html"
    save.plot = TRUE
  }
  
  return(list("object" = object,
              "color.variable" = color.variable,
              "label.variable" = label.variable,
              "filename" = filename,
              "save.plot" = save.plot))
}
