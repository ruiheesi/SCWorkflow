# This code comes from NIDAP 'Dual Labeling [CCBR] [scRNA-seq]' code template

#' @title Plot coexpression of 2 markers using transcript and/or protein expression values 
#' @description Returns individual and combined expression of 2 genes and allows for filtering (optional) of the Seurat object using expression thresholds
#' @details This method provides visualization of coexpression of 2 genes (or proteins) and additional methods for filtering for cells with gene expression values that are above or below thresholds set for one or both markers.  
#' 
#' @param object Seurat-class object
#' @param samples Samples to be included in the analysis 
#' @param marker1 First gene/marker for coexpression analysis
#' @param marker.1.type Slot to use for first marker (choices are "SCT","protein","HTO")
#' @param marker2 Second gene/marker for coexpression analysis
#' @param marker.2.type Slot to use for second marker (choices are "SCT","protein","HTO")
#' @param data.reduction Dimension Reduction method to use for image (options are "umap" or "tsne")
#' @param point.size Point size for image (default is 0.5)
#' @param point.shape Point shape for image (default is 16)
#' @param point.transparency Point transparency for image (default is 0.5)
#' @param add.marker.thresholds Add marker thresholds (default is FALSE)
#' @param marker.1.threshold Threshold set for first marker (default is 0.5)
#' @param marker.2.threshold Threshold set for first marker (default is 0.5)
#' @param filter.data Add new parameter to metadata using marker thresholds (default is FALSE)
#' @param M1.filter.direction Filter to samples that have greater or less than marker 1 threshold (default is "greater than")
#' @param M2.filter.direction Filter to samples that have greater or less than marker 2 threshold (default is "greater than")
#' @param apply.filter.1 If TRUE, apply the first filter (default is TRUE)
#' @param apply.filter.2 If TRUE, apply the second filter (default is TRUE)
#' @param filter.condition If TRUE, apply both filters 1 and 2 and take intersection. If FALSE, apply filters and take the union.
#' @param parameter.name Name for metadata column for new marker filters. (Default is "Marker")
#' @param trim.marker.1 Trim top and bottom percentile of marker 1 signal to pre-scale trim values (below) to remove extremely low and high values (Default is TRUE)
#' @param trim.marker.2 Trim top and bottom percentile of marker 2 signal to pre-scale trim values (below) to remove extremely low and high values (Default is TRUE)
#' @param pre.scale.trim Set trimming percentile values (Defalut is 0.99)
#' @param density.heatmap Creates a additional heatmap showing the density distribution of cells. (Default is FALSE)
#' @param display.unscaled.values Set to TRUE if you want to view the unscaled gene/protein expression values (Default is FALSE)

#' @import Seurat
#' @import ggplot2 
#' @importFrom scales rescale
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid grid.draw
#' @importFrom dplyr arrange mutate case_when
#' @importFrom stats quantile
#'   
#' @export 
#' 
#' @return a seurat object with optional additional metadata for cells that are positive or negative for gene markers and a coexpression plot.

DualLabeling <- function(object,
                         samples,
                         marker1,
                         marker.1.type,
                         marker2,
                         marker.2.type,
                         data.reduction = "umap",
                         point.size = 0.5,
                         point.shape = 16,
                         point.transparency = 0.5,
                         add.marker.thresholds = FALSE,
                         marker.1.threshold = 0.5,
                         marker.2.threshold = 0.5,
                         filter.data = FALSE,
                         M1.filter.direction = "greater than",
                         M2.filter.direction = "greater than",
                         apply.filter.1 = TRUE,
                         apply.filter.2 = TRUE,
                         filter.condition = TRUE,
                         parameter.name = "Marker",
                         trim.marker.1 = TRUE,
                         trim.marker.2 = TRUE,
                         pre.scale.trim = 0.99,
                         density.heatmap = FALSE,    
                         display.unscaled.values = FALSE){

  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##

  if(!(marker1 %in% rownames(object))){
    stop(paste0("ERROR: ",marker1," is not found in dataset"))
  }
  if(!(marker2 %in% rownames(object))){
    stop(paste0("ERROR: ",marker2," is not found in dataset"))
  }
  if(!(marker.1.type %in% names(object@assays))){
    stop(paste0("ERROR: ",marker.1.type," slot is not found in dataset"))
  }
  if(!(marker.2.type %in% names(object@assays))){
    stop(paste0("ERROR: ",marker.2.type," slot is not found in dataset"))
  }
    
  
## --------------- ##
## Functions       ##
## --------------- ##
  
#Function for drawing overlay images for umap/tsne:
gg.overlay <- function(so.sub,df,marker1,marker2){
  df %>% arrange(mark1.scale) -> df    
  xmin = min(df$dr1) - 0.1*min(df$dr1)
  xmax = max(df$dr1) + 0.1*min(df$dr1)
  
  gg.z1 <- ggplot(df, aes(dr1,dr2))+
    geom_point(color=rgb(red=df$mark1.scale,green=0,blue=0),shape=point.shape,size=point.size, alpha=point.transparency)+
    theme_classic() +
    xlab(paste0(data.reduction,"-1")) + 
    ylab(paste0(data.reduction,"-2")) +
    ggtitle(marker1) +
    coord_fixed()
  
  df %>% arrange(mark2.scale) -> df
  
  gg.z2 <- ggplot(df, aes(dr1,dr2))+
    geom_point(color=rgb(red=0,green=df$mark2.scale,blue=0),shape=point.shape,size=point.size, alpha=point.transparency)+
    theme_classic() +
    xlab(paste0(data.reduction,"-1")) + 
    ylab(paste0(data.reduction,"-2")) +
    ggtitle(marker2) +
    coord_fixed()
  
  df %>% mutate(avg = mark2.scale+mark1.scale) %>% arrange(avg) -> df
  
  gg <- ggplot(df, aes(dr1,dr2))+
    geom_point(color=rgb(red=df$mark1.scale,green=df$mark2.scale,blue=0),shape=point.shape,size=point.size, alpha=point.transparency)+
    theme_classic() +
    xlab(paste0(data.reduction,"-1")) + 
    ylab(paste0(data.reduction,"-2")) +
    ggtitle("Combined") +
    coord_fixed()
  
  return(list(gg.z1,gg.z2,gg))
}

#Function for plotting expression data in xy overlay format: 
gg.overlay2 <- function(so.sub,df,marker1,marker2){
  df %>% arrange(mark1.scale) -> df    
  # Create unscaled axis labels
  
  if(display.unscaled.values == TRUE){
    label1_min <- paste("unscaled min:", round(min(mark1),digits = 2))
    label1_max <- paste("unscaled max:", round(max(mark1),digits = 2))
    label1 <- paste(as.character(marker1), label1_min, label1_max, sep = "\n")
    
    label2_min <- paste("unscaled min:", round(min(mark2),digits = 2))
    label2_max <- paste("unscaled max:", round(max(mark2),digits = 2))
    label2 <- paste(as.character(marker2), label2_min, label2_max, sep = "\n")} else {
      
      label1 <- as.character(marker1)
      label2 <- as.character(marker2)
    }
  
  gg.z1 <- ggplot(df, aes(mark1.scale,mark2.scale))+
    geom_point(color=rgb(red=df$mark1.scale,green=0,blue=0),shape = 20,size=point.size)+
    theme_classic() +
    xlab(label1) + 
    ylab(label2) +
    coord_fixed()
  
  df %>% arrange(mark2.scale) -> df
  
  gg.z2 <- ggplot(df, aes(mark1.scale,mark2.scale))+
    geom_point(color=rgb(red=0,green=df$mark2.scale,blue=0),shape = 20,size=point.size)+
    theme_classic() +
    xlab(label1) + 
    ylab(label2) +
    coord_fixed()
  
  df %>% mutate(avg = mark2.scale+mark1.scale) %>% arrange(avg) -> df
  
  gg <- ggplot(df, aes(mark1.scale,mark2.scale))+
    geom_point(color=rgb(red=df$mark1.scale,green=df$mark2.scale,blue=0),shape = 20,size=point.size)+
    theme_classic() +
    xlab(label1) + 
    ylab(label2) +
    coord_fixed()
  
  if(add.marker.thresholds==TRUE){
    gg.z1 <- gg.z1 + 
      geom_vline(xintercept=t1,linetype="dashed") +
      geom_hline(yintercept=t2,linetype="dashed") 
    gg.z2 <- gg.z2 + 
      geom_vline(xintercept=t1,linetype="dashed") +
      geom_hline(yintercept=t2,linetype="dashed") 
    gg <- gg +
      geom_vline(xintercept=t1,linetype="dashed") +
      geom_hline(yintercept=t2,linetype="dashed") 
  }
  
  return(list(gg.z1,gg.z2,gg))
}

## --------------- ##
## Main Code Block ##
## --------------- ##

# Load and subset data using sample names

if("active.ident" %in% slotNames(object)){
  sample_name = as.factor(object@meta.data$orig.ident)
  names(sample_name)=names(object@active.ident)
  object@active.ident <- as.factor(vector())
  object@active.ident <- sample_name
  so.sub = subset(object, ident = samples)
} else {
  sample_name = as.factor(object@meta.data$orig.ident)
  names(sample_name)=names(object@active.ident)
  object@active.ident <- as.factor(vector())
  object@active.ident <- sample_name
  so.sub = subset(object, ident = samples)
}

t1 = marker.1.threshold
t2 = marker.2.threshold

#Select marker 1 values and scale:
mark1 = so.sub@assays[[marker.1.type]]@scale.data[marker1,]
if(trim.marker.1 == TRUE){
  q1 = quantile(mark1,pre.scale.trim)
  q0 = quantile(mark1,1-pre.scale.trim)
  mark1[mark1<q0]=q0
  mark1[mark1>q1]=q1
}
mark1.scale <- rescale(mark1, to=c(0,1))

#Select marker 2 values and scale:
mark2 = so.sub@assays[[marker.2.type]]@scale.data[marker2,]
if(trim.marker.2 == TRUE){
  q1 = quantile(mark2,pre.scale.trim)
  q0 = quantile(mark2,1-pre.scale.trim)
  mark2[mark2<q0]=q0
  mark2[mark2>q1]=q1
}
mark2.scale <- rescale(mark2, to=c(0,1))

#Draw Plots:
df <- data.frame(cbind(dr1=so.sub@reductions[[data.reduction]]@cell.embeddings[,1],
                       dr2=so.sub@reductions[[data.reduction]]@cell.embeddings[,2],
                       mark1.scale,mark2.scale))
gg.list <- gg.overlay(so.sub,df,marker1,marker2)
gg.list2 <- gg.overlay2(so.sub,df,marker1,marker2)

# n=3
# imageWidth = 1000*n
# imageHeight = 1000*2
# dpi = 300
# 
# png(
#   width=imageWidth,
#   height=imageHeight,
#   units="px",
#   pointsize=4,
#   bg="white",
#   res=dpi,
#   type="cairo")

if(density.heatmap == TRUE){
  x <- df$mark1.scale
  y <- df$mark2.scale
  
  df_heatmap <- data.frame(x = x, y = y, 
      d = densCols(x, y, nbin = 1000, 
      colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
  
  p <- ggplot(df_heatmap) +
    geom_point(aes(x, y, col = d), size = 1) +
    scale_color_identity() + xlab(marker1) + ylab(marker2) +
    theme_bw() + geom_vline(xintercept=t1,linetype="dashed") +
    geom_hline(yintercept=t2,linetype="dashed") 
  
  
  grob <- arrangeGrob(gg.list[[1]], gg.list[[2]], gg.list[[3]], gg.list2[[1]],gg.list2[[2]],gg.list2[[3]],p,ncol=3)
} else {
  grob <- arrangeGrob(gg.list[[1]], gg.list[[2]], gg.list[[3]], gg.list2[[1]],gg.list2[[2]],gg.list2[[3]],ncol=3)
}

#Applying Filters to Data using Thresholds:
if(filter.data == TRUE){
  df <- df %>% mutate(sample=so.sub@meta.data$sample) %>% 
    mutate(cellbarcode=rownames(so.sub@meta.data))
  
  if(M1.filter.direction == "greater than"){
    ind1 <- df$`mark1.scale` > t1
  } else {
    ind1 <- df$`mark1.scale` < t1
  }
  
  cat("\n")
  print("Marker 1 filter:")
  print(sum(ind1))
  
  if(M2.filter.direction == "greater than"){
    ind2 <- df$`mark2.scale` > t2
  } else {
    ind2 <- df$`mark2.scale` < t2
  }
  
  cat("\n")
  print("Marker 2 filter:")
  print(sum(ind2))
  
  if(apply.filter.1){
    if(apply.filter.2){
      if(filter.condition == TRUE){
        df <- df[c(ind1&ind2),]
      } else {
        df <- df[c(ind1|ind2),]  
      }  
    } else {
      df <- df[ind1,]   
    }
  } else {
    if(apply.filter.2){
      df <- df[ind2,]    
    }
  }
  
  colnames(df)[3:4]= c(marker1,marker2)
  so.sub@meta.data %>% 
    mutate(x = case_when(rownames(so.sub@meta.data) %in% df$cellbarcode ~ TRUE, TRUE ~ FALSE)) -> so.sub.df
  colnames(so.sub@meta.data) <- sub("x",parameter.name,colnames(so.sub@meta.data))
  
  cat("\n")
  print("Final Breakdown:")
  print(addmargins(table(so.sub.df[[parameter.name]],so.sub.df$sample_name)))
  rownames(so.sub.df) <- rownames(so.sub@meta.data)
  so.sub@meta.data <- so.sub.df
  
  cat("\n")
  print("After Filter applied:")
  print(dim(df)[1])
}

result.list <- list("so" = so.sub,"plot" = grob)

return(result.list)
}