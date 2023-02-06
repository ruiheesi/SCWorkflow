# This code comes from NIDAP 'Harmony Batch Correction [scRNA-seq][CCBR]' code template
# Template Manager https://nidap.nih.gov/workspace/vector/templates/ri.vector.main.template.a097334e-21ac-466c-b467-c761b2c25f02
# Documentation https://nidap.nih.gov/workspace/notepad/view/ri.notepad.main.notepad.facd9ed0-90fc-4c55-bbb8-0153f06858d4

#' @title Performs batch correction on Seurat-class object
#' @description Returns a Seurat-class object with adjusted cell embeddings and gene expression data to account for variations in sample collection
#' @details Takes in a list of genes inputted by the user, displays gene expression information in particular slot-assay with (optional) outliers removed
#' 
#' @param seurat_object Seurat-class object
#' @param variable_features Number of most variable genes to subset the gene expression data by
#' @param genes_to_add Add genes that might not be found among variably expressed genes
#' @param group.by.vars which variable should be accountted for when running batch correction

#' @import Seurat 
#' @import harmony
#' @import gridExtra
#' @import RColorBrewer
#' @import ggplot2
#'   
#' @export

#' @return Seurat-class Object with Harmony-adjusted gene expression (SCT slot) and tSNE cell embedding


harmony_batch_correct <- function(so, 
                                  variable_features = 2000, 
                                  genes_to_add = c(),
                                  group.by.vars) {
  
  ## -------------------------- ##
  ## Error and Warning Messages ##
  ## -------------------------- ##
  
  if(is.null(genes_to_add)){
    print("no genes will be added")
  } else if (all(!genes_to_add %in% rownames(so))){
    warning("specified genes were not found and therefore cannot be added")
  }
  
  if (variable_features > length(VariableFeatures(so))){
    stop("Number of variable features to subset by cannot exceed the total number of variable genes in the data")
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  seur.SCT <- so@assays$SCT@scale.data
  
  genes_of_interest <- genes_to_add[genes_to_add %in% rownames(so@assays$SCT@scale.data)]
  
  # Add more genes to analyze
  VariableFeatures(so) <- c(VariableFeatures(so), genes_of_interest) 
  mvf <- VariableFeatures(so)[1:(variable_features + length(genes_of_interest))] 
  
  # Check if variable genes are found in scale data genes
  mvf <- mvf[mvf %in% rownames(seur.SCT)]
  
  print(dim(so@assays$SCT@scale.data))
  seur.SCT <- seur.SCT[mvf,]
  seur.SCT <- t(seur.SCT) 
  # need the original loadings and embeddings to compare with my SVD
  seur.loads <- so@reductions$pca@feature.loadings
  seur.pca <- so@reductions$pca@cell.embeddings
  sm.pca <- seur.pca[1:10, 1:10]
  sm.seur.ldngs <- seur.loads[1:10, 1:10]
  
  # SVD on scaled counts
  
  variable_features <- length(mvf)
  Sys.time()
  pppca <- svd(seur.SCT) 
  Sys.time() 
  ppembed <- pppca$u %*% diag(pppca$d)
  pcnames <- vector(mode = "character")
  for (i in 1:dim(ppembed)[2])pcnames[i] <- paste("PC", i, sep = "_")
  colnames(ppembed) <- pcnames
  rownames(ppembed) <- rownames(seur.SCT)
  sm.ppembed <- ppembed[1:10, 1:10]
  ppldngs <- pppca$v
  colnames(ppldngs) <- pcnames
  rownames(ppldngs) <- mvf 
  sm.ppldngs <- ppldngs[1:10, 1:10]
  ppembed[1:10, 1:10]
  ppldngs[1:10, 1:10]
  
  so@reductions$pca@cell.embeddings <- ppembed
  so@reductions$pca@feature.loadings <- ppldngs
  so@reductions$pca@stdev <- pppca$d
  print(Sys.time())
  
  so <- RunHarmony(so, group.by.vars,
                   do_pca=FALSE,
                   assay.use = "SCT",
                   plot_convergence = TRUE,
                   return_object=TRUE)
  
  head(so@reductions$harmony@cell.embeddings)
  head(so@reductions$harmony@feature.loadings)
  
  so <- RunUMAP(so, reduction = "harmony",dims=1:4)
  so <- RunTSNE(so, reduction = "harmony",dims=1:4)
  
  sdat <- data.frame(as.vector(so@reductions$tsne@cell.embeddings[,1]),
                     as.vector(so@reductions$tsne@cell.embeddings[,2]),
                     so@meta.data[eval(parse(text = "group.by.vars"))])
  names(sdat) <- c("TSNE1","TSNE2","ident")
  
  n <- 2e3
  set.seed(10)
  ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
  ourColorSpace <- as(ourColorSpace, "LAB")
  
  
  distinctColorPalette <-function(k=1) {
    currentColorSpace <- ourColorSpace@coords
    # Set iter.max to 20 to avoid convergence warnings.
    set.seed(1)
    km <- kmeans(currentColorSpace, k, iter.max=20)
    colors <- unname(hex(LAB(km$centers)))
    return(colors)
  }   
  
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
  cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  g1 <- ggplot(sdat, aes(x=TSNE1, y=TSNE2)) +
    theme_bw() +
    theme(legend.title=element_blank()) +
    geom_point(aes(colour=ident),size=1) +
    scale_color_manual(values=cols) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position="top",
          panel.background = element_blank(),
          legend.text=element_text(size=rel(0.5))) +
    guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) 
  
  solist <- list(so, ppldngs)
  
  # Calculate adjusted gene expression from embeddings
  so <- solist[[1]]
  ppldngs <- solist[[2]]
  seur.SCT <- t(so@assays$SCT@scale.data)
  harm.embeds <- so@reductions$harmony@cell.embeddings
  #print(dim(harm.embeds))
  harm.lvl.backcalc <- harm.embeds %*% t(ppldngs)
  
  sm.harmexpr <- harm.lvl.backcalc[1:10, 1:4]
  sm.comp.scld <- (seur.SCT[rownames(sm.harmexpr), colnames(sm.harmexpr)])

  # Convert to ggplot that includes color by samples for higher resolution of batch correction, use for loop to screen for different genes
  so@assays$SCT@scale.data <- t(harm.lvl.backcalc)

  return(so)
}
