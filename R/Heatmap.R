# This code comes from NIDAP 'Heatmap for Single Cell Data [scRNA-Seq][CCBR]' code template

#' @title Plot coexpression of 2 markers using transcript and/or protein expression values 
#' @description Returns individual and combined expression of 2 genes and allows for filtering (optional) of the Seurat object using expression thresholds
#' @details This method provides visualization of coexpression of 2 genes (or proteins) and additional methods for filtering for cells with gene expression values that are above or below thresholds set for one or both markers.  
#' 
#' @param object Seurat-class object
#' @param sample.names Sample names
#' @param metadata Metadata column to plot
#' @param transcripts Transcripts to plot
#' @param proteins Proteins to plot (default is NULL)
#' @param plot.title Title of plot (default is "Heatmap")
#' @param add.gene.or.protein Add Gene or protein annotations (default is FALSE)
#' @param protein.annotations Protein annotations to add (defulat is NULL)
#' @param rna.annotations Gene annotations to add (default is NULL)
#' @param arrange.by.metadata Arrange by metadata (default is TRUE)
#' @param add.row.names Add row names (default is FALSE)
#' @param add.column.names Add column names (default is FALSE)
#' @param set.seed Seed for colors (default is 6)
#' @param scale.data Perform z-scaling on rows (default is TRUE)
#' @param trim.outliers Remove outlier data (default is TRUE)
#' @param trim.outliers.percentage Set outlier percentage (default is 0.01)
#' @param order.heatmap.rows Order heatmap rows (default is FALSE)
#' @param row.order Set row order to gene vector (default is NULL)
#' 
#' @import Seurat
#' @import ggplot2 
#' @import pheatmap
#' @import dendsort
#' @import dendextend
#' @import colorspace
#'   
#' @export
#'   
Heatmap <- function(object,
                    sample.names,
                    metadata,
                    transcripts,
                    proteins = NULL,
                    plot.title = "Heatmap",
                    add.gene.or.protein = TRUE,
                    protein.annotations = NULL,
                    rna.annotations = NULL,
                    arrange.by.metadata = TRUE,
                    add.row.names = TRUE,
                    add.column.names = FALSE,
                    set.seed = 6,
                    scale.data = TRUE,
                    trim.outliers = TRUE,
                    trim.outliers.percentage = 0.01,
                    order.heatmap.rows = FALSE,
                    row.order = c()) {
  ## -------------------------------- ##
  ## Functions                        ##
  ## -------------------------------- ##
  
  n <- 2e3
  
  set.seed(set.seed)
  ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
  ourColorSpace <- as(ourColorSpace, "LAB")
  
  
  distinctColorPalette <-function(k=1,seed) {
    currentColorSpace <- ourColorSpace@coords
    # Set iter.max to 20 to avoid convergence warnings.
    set.seed(seed)
    km <- kmeans(currentColorSpace, k, iter.max=20)
    colors <- unname(hex(LAB(km$centers)))
    return(colors)
  }   
  
  pal = function (n, h=c(237, 43), c=100, l=c(70, 90), power=1, fixup=TRUE, gamma=NULL, alpha=1, ...) {
    if (n < 1L) 
      return(character(0L))
    h <- rep(h, length.out = 2L)
    c <- c[1L]
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, -1, length = n)
    rval <- hex(
      polarLUV(
        L = l[2L] - diff(l) * abs(rval)^power[2L], 
        C = c * abs(rval)^power[1L],
        H = ifelse(rval > 0, h[1L], h[2L])
      ),
      fixup=fixup, ...
    )
    if (!missing(alpha)) {
      alpha <- pmax(pmin(alpha, 1), 0)
      alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                      width = 2L, upper.case = TRUE)
      rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
  }
  
  #Color selections for heatmap:
  np0 = pal(100)
  np1 = diverge_hcl(100, c=100, l=c(30, 80), power=1)  #Blue to Red
  np2 = heat_hcl(100, c=c(80, 30), l=c(30, 90), power=c(1/5, 2))  #Red to Vanilla
  np3 = rev(heat_hcl(100, h=c(0, -100), c=c(40, 80), l=c(75, 40), power=1)) #Violet to Pink
  np4 = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100))
  np5 = colorRampPalette(c("steelblue","white", "red"))(100) #Steelblue to White to Red
  
  np = list(np0, np1, np2, np3, np4, np5)
  names(np) = c("Default","Blue to Red","Red to Vanilla","Violet to Pink","Bu Yl Rd","Bu Wt Rd")
  
  
  doheatmap <- function(dat, clus, clus2, rn, cn, col) {
    require(pheatmap)
    require(dendsort)
    if (scale.data == TRUE) {
      tmean.scale = t(scale(t(dat)))
      tmean.scale = tmean.scale[is.finite(rowSums(tmean.scale)),]
    } else {
      tmean.scale = dat
    }
    if(trim.outliers == TRUE){
      quantperc <- trim.outliers.percentage
      upperquant <- 1-quantperc
      for(i in 1:nrow(tmean.scale)){
        data <- tmean.scale[i,] 
        dat.quant <- quantile(data,probs=c(quantperc,upperquant))
        data[data > dat.quant[2]] <- dat.quant[2]
        data[data < dat.quant[1]] <- dat.quant[1]
        tmean.scale[i,] <- data    
      }
    }
    col.pal <- np[[col]]
    drows1 <- "euclidean"
    dcols1 <- "euclidean"
    minx = min(tmean.scale)
    maxx = max(tmean.scale)
    treeheight <- 25
    breaks = seq(minx, maxx, length=100)
    legbreaks = seq(minx, maxx, length=5)
    breaks = sapply(breaks, signif, 4)
    legbreaks = sapply(legbreaks, signif, 4)
    
    #Run cluster method: 
    hc = hclust(dist(t(tmean.scale)), method="complete")
    hcrow = hclust(dist(tmean.scale), method="complete")
    
    if (clus) {
      sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))
    } else {
      sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
    }
    
    if (clus2) {
      rowclus <- sort_hclust(hcrow)
    } else {
      rowclus = FALSE
    }
    
    pathname <- stringr::str_replace_all(plot.title, "_", " ") 
    pathname <- stringr::str_wrap(pathname,50)
    
    hm.parameters <- list(
      tmean.scale, 
      color=col.pal,
      legend_breaks=legbreaks,
      scale="none",
      treeheight_col=treeheight,
      treeheight_row=treeheight,
      kmeans_k=NA,
      breaks=breaks,
      fontsize_row=5,
      fontsize_col=10,
      show_rownames=rn, 
      show_colnames=cn,
      main=pathname,
      clustering_method="complete",
      cluster_rows=rowclus, 
      cluster_cols=clus,
      cutree_rows=1,
      clustering_distance_rows=drows1, 
      clustering_distance_cols=dcols1,
      annotation_col = annotation_col,
      annotation_colors = annot_col,
      labels_col = labels_col
    )
    
    mat = t(tmean.scale)
    
    callback = function(hc, mat) {
      dend=rev(dendsort(as.dendrogram(hc)))
      dend %>% dendextend::rotate(c(1:length(dend))) -> dend
      as.hclust(dend)
    }
    do.call("pheatmap", c(hm.parameters, list(clustering_callback=callback)))
  }
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  # load data
  samples_to_include = sample.names
  samples_to_include <- samples_to_include[samples_to_include != ""]
  samples_to_include <- gsub("-","_",samples_to_include)
  
  #Clean up transcript names and print missing genes:
  transcripts = gsub(" ","",transcripts) 
  if(transcripts[1] != ""){
    genesmiss = setdiff(transcripts,rownames(object$SCT@scale.data))
    if(length(genesmiss)>0){
      print(paste("missing genes:", genesmiss))
    }
  }
  transcripts = transcripts[transcripts %in% rownames(object$SCT@scale.data)]
  
  #Clean up protein names and print missing proteins:
  if(!is.null(object@assays$Protein)){
    proteins = gsub(" ","",proteins) 
    if(proteins[1] != ""){
      protmiss = setdiff(proteins,rownames(object$Protein@scale.data))
      if(length(protmiss)>0){
        print(paste("missing proteins:", protmiss))
      }
    }
    proteins = proteins[proteins %in% rownames(object$Protein@scale.data)]
  }
  
  #collect transcript expression data from SCT slot
  df.mat1 = NULL
  if(length(transcripts)>0){
    if(length(transcripts)==1){
      df.mat1 <- vector(mode="numeric",length=length(object$SCT@scale.data[transcripts,]))
      df.mat1 <- object$SCT@scale.data[transcripts,]
    } else {
      df.mat1 <- as.matrix(object$SCT@scale.data[transcripts,])
    }
  }
  
  #collect protein expression values
  df.mat2 = NULL
  if(!is.null(object@assays$Protein)){
    if(length(proteins)>0){
      if(length(proteins)==1){
        df.mat2 <- vector(mode="numeric",length=length(object$Protein@scale.data[proteins,]))
        df.mat2 <- object$Protein@scale.data[proteins,]
        protname <- paste0(proteins,"_Prot")
      } else {
        df.mat2 <- as.matrix(object$Protein@scale.data[proteins,])
        protname <- paste0(proteins,"_Prot")
        rownames(df.mat2) <- protname
      }
    }
  }
  
  #Put together transcript and protein data
  df.mat <- rbind(df.mat1,df.mat2)
  if(!is.null(df.mat1)){
    rownames(df.mat)[rownames(df.mat)=="df.mat1"] <- transcripts
  }
  
  if(!is.null(df.mat2)){
    rownames(df.mat)[rownames(df.mat)=="df.mat2"] <- protname
  }
  
  df.mat <- df.mat[sort(rownames(df.mat)),]
  
  #Set row order to specific genes/proteins  
  if(order.heatmap.rows == TRUE){
    row.order = row.order[row.order %in% rownames(df.mat)]
    row.order = c(row.order,setdiff(rownames(df.mat),row.order))
    df.mat <- df.mat[row.order,]
    clusrows <- FALSE
  } else {
    clusrows <- TRUE
  }
  
  metadata <- sub("orig_ident","orig.ident",metadata)
  #Set up annotation track: 
  object@meta.data %>% dplyr::filter(orig.ident %in% samples_to_include) -> metadata_table
  
  if(!"Barcode" %in% metadata){
    metadata = c(metadata,"Barcode")
  }
  
  metadata_table %>% select(metadata) -> annot
  a = dim(annot)[2] - 1
  
  #Addition of gene and/or protein expression data to annotation tracks:
  if(add.gene.or.protein == TRUE){
    if(length(protein.annotations) > 0){
      annot1 <- as.matrix(object$Protein@scale.data[protein_annotations,])
      if(length(protein.annotations)==1){
        annot1 <- annot1[match(annot$Barcode,rownames(annot1))]
        protname <- paste0(protein.annotations,"_Prot")
      } else {
        annot1 <- annot1[,match(annot$Barcode,colnames(annot1))]
        annot1 <- t(annot1)
        colnames(annot1) = paste0(colnames(annot1),"_Prot")
      }
    }
    
    if(length(rna.annotations)>0){
      annot2 <- as.matrix(object$SCT@scale.data[rna.annotations,])
      if(length(rna.annotations)==1){
        annot2 <- annot2[match(annot$Barcode,rownames(annot2))]
      } else {    
        annot2 <- annot2[,match(annot$Barcode,colnames(annot2))]
        annot2 <- t(annot2)
      }
    }
  }
  
  if(exists("annot1")){
    annot <- cbind(annot,annot1)
    colnames(annot)[colnames(annot)=="annot1"] <- protname
  }
  if(exists("annot2")){
    annot <- cbind(annot,annot2)
    colnames(annot)[colnames(annot)=="annot2"] <- rna.annotations
  }
  
  #Arrange columns by metadata tracks: 
  if(arrange.by.metadata == TRUE){
    annot %>% arrange(across(all_of(metadata))) -> annot
    df.mat <- df.mat[,match(annot$Barcode,colnames(df.mat))] 
    df.mat <- df.mat[ , apply(df.mat, 2, function(x) !any(is.na(x)))]
    cluscol = FALSE
  } else {
    cluscol = TRUE
  }
  
  #Set up annotation track and colors:
  annotation_col = as.data.frame(unclass(annot[,!(names(annot) %in% "Barcode")]))
  annotation_col %>% mutate_if(is.logical, as.factor) -> annotation_col
  rownames(annotation_col) <- annot$Barcode
  if(dim(annot)[2] == 2){
    annottitle = colnames(annot)[1]
    colnames(annotation_col) = annottitle
  }
  annot_col = list()
  groups=colnames(annotation_col)
  
 q=0
  for(i in 1:dim(annotation_col)[2]){
    q = q+length(unique(annotation_col[,i]))
  } 
  
  colors=distinctColorPalette(q,5)
  
  b=1
  i=1
  nam = NULL
  col <- NULL
  annot_col <- NULL
  for (i in 1:length(groups)){
    nam <- groups[i]
    if(class(annotation_col[,i]) != "numeric"){
      grp <- as.factor(annotation_col[,i])
      c <- b+length(levels(grp))-1
      col = colors[b:c]
      names(col) <- levels(grp)
      assign(nam,col)
      annot_col = append(annot_col,mget(nam))
      b = c+1
      i=i+1
    }
    else{
      grp <- annotation_col[,i]
      np5 = colorRampPalette(c("steelblue","white", "red"))(length(grp))
      col=np5
      assign(nam,col)
      annot_col = append(annot_col,mget(nam))
    }
  }
  
  print(paste0("The total number of genes in heatmap: ", nrow(df.mat)))
  labels_col <- colnames(df.mat)
  
  #Draw heatmap
  p = doheatmap(dat=df.mat, clus=cluscol, clus2=clusrows, rn=add.row.names, cn=add.column.names, col="Bu Yl Rd")
  
  #Return expression matrix used in heatmap
  as.data.frame(df.mat) %>% rownames_to_column("gene") -> heatmap.df
  heatmap.res <- list("plot" = p, "data" = heatmap.df)
  return(heatmap.res)
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.


 