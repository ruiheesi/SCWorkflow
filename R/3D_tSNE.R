# This code comes from NIDAP 'Plot 3D TSNE from Seurat Object [scRNA-Seq][CCBR]' code template

#' @title Plot 3D TSNE from Seurat Object and saves as html file
#' @description Returns 3D tSNE plot for seurat object
#' @details This method provides visualization of 3D tSNE plot
#' 
#' @param object Seurat-class object
#' @param color.variable Metadata column to use for color
#' @param label.variable Metadata column to use for label
#' @param fileName Filename for saving plot (default is "plot.html")
#' @param npcs Number of PC's (default is 15)
#' @param save.plot Save plot (default is FALSE)
#' 
#' @import Seurat
#' @import ggplot2 
#' @importFrom plotly as_widget plot_ly
#'   
#' @export

tSNE3D <- function(object,
                   color.variable,
                   label.variable,
                   fileName = "plot.html",
                   save.plot = FALSE,
                   npcs = 15){
  
  color.variable <- sub("orig_ident","orig.ident",color.variable)
  label.variable <- sub("orig_ident","orig.ident",label.variable)
  cols=c("blue","green","red","orange","purple4","darkcyan","magenta2",
         "darkred","darkorange")
  
  #Run TSNE again to get 3d coordinates:
  object <- RunTSNE(object, assay="SCT",
                dims=1:npcs,
                dim.embed = 3, 
                seed.use=1)
  
  tsne.coord <- as.data.frame(object@reductions$tsne@cell.embeddings)
  
  tsne.df <- data.frame(TSNE1=tsne.coord$tSNE_1,
                        TSNE2=tsne.coord$tSNE_2,
                        TSNE3=tsne.coord$tSNE_3,
                        colorvar=as.factor(object@meta.data[[color.variable]]),
                        label = object@meta.data[[label.variable]])
  
  fig <- plot_ly(tsne.df, x = ~TSNE1, 
                 y = ~TSNE2, 
                 z = ~TSNE3, 
                 color = ~colorvar, 
                 colors = cols,
                 type="scatter3d",
                 mode="markers",
                 hoverinfo = 'text',
                 text = ~label,
                 size= 4) 
  
  legend=TRUE
  if(legend==FALSE){
    fig <- hide_legend(fig)
  }
  
  if(save.plot == TRUE){
    htmlwidgets::saveWidget(as_widget(fig), fileName, selfcontained=TRUE)
  }
  
  tsne.results <- list("plot"  = fig, "data" = tsne.df)
  return(tsne.results) 
}
