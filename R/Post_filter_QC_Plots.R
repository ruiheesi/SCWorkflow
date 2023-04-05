#' @title Post-filter QC Plots
#' @description This set of plots compares several QC parameters across samples 
#' after initial QC and filtering. 
#' @details This is Step 2 in the basic Single-Cell RNA-seq workflow. 
#' This Function is used to evaluate different factors effecting data quality 
#' that should be taken into consideration for downstream processing
#' @param object Please input a filtered Seurat Object. 
#' This should be the output from the Initial QC template
#' @param image.type Remember that svgs are much larger than pngs, 
#' so we recommend doing everything first in png, then rerunning to output 
#' specific svgs as needed
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
#' @export
#' 
#' @return Original Seurat Object and QC plots
postFilterQC <- function(object,
                         image.type = 'png'
){
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  .qFilter <- function(x){	
    qc.df %>% dplyr::filter(variable == x)
  }
  
  ## Create Histogram Plot
  .plotHist <- function(count.df){	
    g=ggplot(count.df) + 
      theme_bw() +
      geom_density(aes(x = value, colour = orig.ident)) +
      labs(x = NULL) +
      theme(legend.position='right',legend.text=element_text(size=10),
            legend.title=element_blank()) + 
      ggtitle(count.df$variable[1]) +
      scale_x_continuous(trans='log10') + 
      scale_color_manual(values = col3) 
    return(g)
  }
  
  ## Create Violin Plot
  .plotViolin <- function(count.df){
    axis.lab = unique(count.df$orig.ident)	
    
    v<-ggplot(count.df, aes(x=orig.ident, y=value)) +
      ggtitle(count.df$variable[1]) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.text=element_text(size=rel(1)),
            legend.title=element_blank(), 
            axis.text=element_text(size=10),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            #axis.text.x=element_text(angle=45,hjust=1),
            plot.title = element_text(size = 20, face = "bold")) +
      geom_violin(aes(fill=as.factor(orig.ident))) +  
      scale_fill_manual(values = col3) +
      geom_boxplot(width=.1) +
      scale_x_discrete(limits = as.vector(axis.lab)) 
    return(v)
    
    
  }
  
  ## Create Scatter Plot
  .plotScatter <- function(count.df,counts){	
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
    return(g)
  }
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  col1<-brewer.pal(8, "Set3")[-2] 
  col2<-c(col1,brewer.pal(8,"Set2")[3:6])
  col3<-c(col2,"#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","
          #f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324","#800000",
          "#aaffc3","#808000","#000075","#a9a9a9","#808080","#A9A9A9","#8B7355") 
  
  
  
  
  #in case you want to this on a merged object
  if (class(object) =="Seurat") {
    x =list()
    x[[1]] <- object
    object <- x
  }
  
  ## Metadata Info
  all.columns <- unique(unlist(
    sapply(
      seq_along(object), 
      function(i) colnames(object[[i]]@meta.data))))
  
  qc.df <- array(0,dim=c(0,3))
  for (i in 1:length(object)){
    so <- object[[i]]
    
    #Add missing columns to metadata
    missing.columns <- setdiff(all.columns,colnames(so@meta.data))
    
    for(i in missing.columns){
      so <- AddMetaData(so,rep(0,ncol(so)), i)
    }
    df.m <- melt(so@meta.data)%>%suppressMessages()
    qc.df <- rbind(qc.df,df.m)
  }
  
  ###################################
  ## Create Images
  
  qc.count <- lapply(unique(qc.df$variable), function(x) {.qFilter(x)})
  qc.count[[1]] %>% 
    dplyr::filter(variable=="nCount_RNA") %>% 
    pull(value) -> RNAcounts
  
  grobs <- lapply(seq_along(qc.count), 
                  function(x){arrangeGrob(grobs = list(
                    .plotScatter(qc.count[[x]],RNAcounts),
                    .plotHist(qc.count[[x]]),
                    .plotViolin(qc.count[[x]])),
                    nrow=1,ncol=3)})%>%suppressWarnings()
  
  #############################
  ## Plot Image 
  
  grobs <- arrangeGrob(grobs=grobs,nrow=length(grobs))
  
  
  return(list(object=object,plot=grobs))
  
  
  
}