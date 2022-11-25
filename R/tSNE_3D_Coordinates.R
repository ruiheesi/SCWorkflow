#' 3D tSNE Coordinate 
#' Template from v 75
#' 
#' 
#' @param Combine_and_Renormalize Seurat Object output from previous combine and renormalize template/node. 
#' @param Gene Default values: Gapdh, GAPDH. Please enter genes which you would like to visualize. If you don't know what genes to look at, consider consulting your DEG table for interest contrasts. 
#' @param Assay Default value: SCT. Select from SCT,RNA,integrated. Select Assay to Plot (default is SCT).
#' @param Max_sample Default value: 10000. Random subsampling of cells without replacement.  At the moment only ten thousand cells can be displayed in 3d tSNE viewer.
#' @param On_NIDAP Default value: false. This option is to use the API calls
#' 
#' @import Seurat 
#' @import httr
#' @import jsonlite 
#' @import plyr 
#' 
#' @export
#' 
#' @return return the  Register the output dataset on phonograph2 and ret



tSNE_3D_Coordinates <- function(Combine_and_Renormalize,
                                Gene = c("Gapdh","GAPDH"),
                                Assay = "SCT",
                                Max_sample = 10000,
                                On_NIDAP = FALSE) {
  
  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  
  #Basic Parameters
  so <- Combine_and_Renormalize
  
  gene = Gene 
  max_sample<- Max_sample
  
  if (Assay == "SCT"){
    Assay_selected <- so$SCT
      
  }else if (Assay == "RNA"){
    Assay_selected <- so$RNA
    
  }else if ("Assay == integrated"){
    Assay_selected <- so$integrated
    
  }else{
    stop("Invalid Assay value")
  }
  

  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##  
  # Data processing script
  # Select missing gene from dataset
  nogene = gene[!gene %in% rownames(Assay_selected@scale.data)]
  
  if(!is.null(nogene)){
    print("Gene(s) missing from dataset:")
    print(nogene)
  }
  
  # Extract data for selected genes 
  gene = gene[gene %in% rownames(Assay_selected@scale.data)]
  
  gene.mat=Assay_selected@scale.data[gene,]
  
  gene.quant=quantile(gene.mat[gene.mat>1],probs=c(.1,.5,.90))
  gene.mat[gene.mat>gene.quant[3]]=gene.quant[3]
  gene.mat[gene.mat<gene.quant[1]]=0
  
  gene.mat<-as.data.frame(gene.mat)
  colnames(gene.mat)<-gene
  
  # Run TSNE on the imported dataset
  yourseuratobject <- RunTSNE(object = so, reduction = "pca", dim.embed = 3, seed.use = 1)
  
  # Gather spatial coordinates for cells
  x <- yourseuratobject[["tsne"]]@cell.embeddings[,1]
  y <- yourseuratobject[["tsne"]]@cell.embeddings[,2]
  z <- yourseuratobject[["tsne"]]@cell.embeddings[,3]
  
  # Generate raw output dataset with coordinates assinged to selected genes 
  pk=1:length(x) 
  
  if(length(gene)>0){
    rv <- as.data.frame(cbind(pk=pk,x=x,y=y,z=z,so@meta.data,gene.mat))
  }else{
    rv <- as.data.frame(cbind(pk=pk,x=x,y=y,z=z,so@meta.data))
  }
  
  # Checking if the raw output dataset has length does not exceed the user defined max sample number   
  if(nrow(rv) < max_sample){
    ms = nrow(rv)
    
  }else{
    ms = max_sample
  }
  
  #Check for boolean datatype and conver to strings
  for (col_num in 1:ncol(rv)) {        
    if (class(rv[,col_num]) == "logical") {
      rv[,col_num] <- as.character(rv[,col_num])
    }
  }
  
  # rvs is the final output dataset for plotting 
  rvs=rv[sample(nrow(rv),ms),]
  
  return(as.data.frame(rvs))
}
  
