getParam3D <- function(data) {
  if (data == "TEC") {
    object <- selectSRObject("TEC")
    color.variable = "orig.ident"
    label.variable = "immgen_main"
    filename = "output/TEC_tsneplot.html"
    save.plot = TRUE
  } else if (data == "Chariou") {
    object <- selectSRObject("Chariou")
    color.variable = "mouseRNAseq"
    label.variable = "Barcode"
    filename = "output/Chariou_tsneplot.html"
    save.plot = TRUE
  } else if (data == "pbmc-single") {
    object <- selectSRObject("pbmc-single")
    color.variable = "SCT_snn_res.0.2"
    label.variable = "BP_encode_main"
    filename = "output/pbmc-single_tsneplot.html"
    save.plot = TRUE
  } else if (data == "nsclc-multi") {
    object <- selectSRObject("nsclc-multi")
    color.variable = "orig.ident"
    label.variable = "BP_encode_main"
    filename = "output/nsclc-multi_tsneplot.html"
    save.plot = TRUE
  } else if (data == "BRCA") {
    object <- select_crobject("BRCA")
    color.variable = "orig.ident"
    label.variable = "SCT_snn_res.0.4"
    filename = "output/BRCA_tsneplot.html"
    save.plot = TRUE
  }
  
  return(
    list(
      "object" = object,
      "color.variable" = color.variable,
      "label.variable" = label.variable,
      "filename" = filename,
      "save.plot" = save.plot
    )
  )
}
