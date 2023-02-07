#' Post-filter QC Plots template in single-cell-rna-seq-r4 NIDAP environment
#' from v47
#' 
#' @title Post-filter QC Plots
#' @description This set of plots compares several QC parameters across samples after initial QC and filtering. This template is Step 2 in the basic Single-Cell RNA-seq workflow. 
#' @details 
#'
#' @param Seurat_Object Please input a filtered Seurat Object. This should be the output from the Initial QC template.
#' @param Parallelize_Computation Toggle true to parallelize your computations using Spark. Recommended for large datasets.
#' @param Image_type Remember that svgs are much larger than pngs, so we recommend doing everything first in png, then rerunning to output specific svgs as needed.
#' 
#' 
#' @import Seurat
#' @import tidyverse
#' @import ggplot2
#' @import gtable
#' @import gridExtra 
#' @import reshape2
#' @import RColorBrewer
#' @import dplyr
#' @import svglite
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom digest digest
#' @importFrom svglite svglite
#' 
#' 
#' @export
#' 
#' @return Seurat Objects and QC plots that compares several QC parameters across samples after initial QC and filtering. This template is Step 2 in the basic Single-Cell RNA-seq workflow.
  
<<<<<<< HEAD

=======
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
Post_filter_QC <- function(Seurat_Object,
                           Parallelize_Computation = F,
                           Image_type = 'png'
                           ){
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  qfilter <- function(x){
    qc.df %>% dplyr::filter(variable == x)
  }
  
  plothist <- function(count.df){
    g=ggplot(count.df) + 
      theme_bw() +
      geom_density(aes(x = value, colour = orig.ident)) +
      labs(x = NULL) +
      theme(legend.position='right',legend.text=element_text(size=10),
            legend.title=element_blank()) + 
      ggtitle(count.df$variable[1]) +
      scale_x_continuous(trans='log10') + 
      scale_color_manual(values = col3) 
    #scale_linetype_manual(values=rep(c('solid', 'dashed','dotted'),6))
    return(g)
    #print(g)
    #return(NULL)
  }
  
  plotviolin <- function(count.df){
    axislab = unique(count.df$orig.ident)
    
    v=ggplot(count.df, aes(x=orig.ident, y=value)) +
      ggtitle(count.df$variable[1]) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),legend.text=element_text(size=rel(1)),
            legend.title=element_blank(), axis.text=element_text(size=10),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            #axis.text.x=element_text(angle=45,hjust=1),
            plot.title = element_text(size = 20, face = "bold")) +
      geom_violin(aes(fill=as.factor(orig.ident))) +  
      scale_fill_manual(values = col3) +
      geom_boxplot(width=.1) +
      #labs(colour = n, y=m) +
      #geom_jitter(height = 0, width = 0.1, size = 0.1) +
      #scale_y_continuous(trans='log2') + 
      scale_x_discrete(limits = as.vector(axislab)) 
    return(v)
    #print(v)
    #return(NULL)
  }
  
  plotscatter <- function(count.df,counts){
    count.df %>% dplyr::mutate("value2"=counts) -> count.df 
    ylab = as.character(unique(count.df$variable))
    xlab = "RNA Count"
    name = paste(ylab,"vs.",xlab)          
    g <- ggplot(count.df, aes(x=value2, y=value,color = orig.ident)) +
      geom_point(size = 0.5) + 
      theme_classic() +
      theme(legend.position='right',legend.text=element_text(size=10),
            legend.title=element_blank()) + 
      guides(colour = guide_legend(override.aes = list(size=2))) +
      scale_color_manual(values = col3) +
      labs(title=name, x = xlab, y = ylab)
    ggtitle(name) 
    #print(g)
    #return(NULL)
    return(g)
  }
  
  
  col1=brewer.pal(8, "Set3")[-2] 
  col2=c(col1,brewer.pal(8,"Set2")[3:6])
  col3=c(col2,"#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324","#800000","#aaffc3","#808000","#000075","#a9a9a9","#808080","#A9A9A9","#8B7355")
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  ###################################
  ## Create SO object
  
  # load data
  # object.class <- getClass(class(Seurat_Object))
  # 
  # if(object.class@className == "RFoundryObject") {
  #   cat("1. Reading Seurat Object from dataset: RObjectdata.rds\n\n")
  #   
<<<<<<< HEAD

=======
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  #   SO = Seurat_Object$value
  # } else {
  #   cat("1. Reading Seurat Object from dataset: seurat_object.rds\n\n")
  #   
  #   fs <- Seurat_Object$fileSystem()
  #   path <- fs$get_path("seurat_object.rds", 'r')
<<<<<<< HEAD

=======
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  #   SO <- readRDS(path)
  #   
  # }
  
<<<<<<< HEAD
=
  SO <- Seurat_Object

=======
  SO <- Seurat_Object
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  
  #in case you want to redo this on a merged SO
  if (class(SO) =="Seurat") {
    x =list()
    x[[1]] <- SO
    SO <- x
  }
  
  all.columns <- unique(unlist(sapply(seq_along(SO), function(i) colnames(SO[[i]]@meta.data))))
  qc.df <- array(0,dim=c(0,3))
  for (i in 1:length(SO)){
    so <- SO[[i]]
    #print(so)
    
    #Add missing columns to metadata
    missing.columns <- setdiff(all.columns,colnames(so@meta.data))
    
    for(i in missing.columns){
      so <- AddMetaData(so,rep(0,ncol(so)), i)
    }
    df.m <- melt(so@meta.data)%>%suppressMessages()
    qc.df <- rbind(qc.df,df.m)
  }
  
  ###################################
  ## Create Immages
  
  # if (Parallelize_Computation) {
  #   qc.count <- spark.lapply(unique(qc.df$variable), function(x) {qfilter(x)})
  #   qc.count[[1]] %>% dplyr::filter(variable=="nCount_RNA") %>% pull(value) -> RNAcounts
  # } else {
    qc.count <- lapply(unique(qc.df$variable), function(x) {qfilter(x)})
    qc.count[[1]] %>% dplyr::filter(variable=="nCount_RNA") %>% pull(value) -> RNAcounts
  # }
  grobs <- lapply(seq_along(qc.count), function(x) {arrangeGrob(grobs = list(plotscatter(qc.count[[x]],RNAcounts),plothist(qc.count[[x]]),plotviolin(qc.count[[x]])),nrow=1,ncol=3)})%>%suppressWarnings()
  
  #############################
  ## Plot Image 
  
<<<<<<< HEAD

=======
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  ## Set Image Size   
  imageWidth = 5000
  imageHeight = 1000*length(grobs)
  dpi = 300
  
  if (Image_type == 'png') {
    png(
      # filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
  } else {
    svglite::svglite(
      # file=graphicsFile,
      width=round(imageWidth/dpi,digits=2),
      height=round(imageHeight/dpi,digits=2),
      pointsize=1,
      bg="white")
  }
<<<<<<< HEAD

=======
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  
  
  # grobs= grid.arrange(grobs = grobs, nrow = length(grobs))
  grobs = arrangeGrob(grobs=grobs,nrow=length(grobs))
  
  
  # so@meta.data %>% rownames_to_column("Barcode") -> meta.df
  
<<<<<<< HEAD

  cat("\nReturn objects checksum:\n")
  print(digest::digest(so))

=======
  cat("\nReturn objects checksum:\n")
  print(digest::digest(so))
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  
  return(list(so=SO,plot=grobs))
         
  # output <- new.output()
  # output_fs <- output$fileSystem()
  # saveRDS(SO, output_fs$get_path("seurat_object.rds", 'w'))
  # 
  # return(output_fs)
  
}