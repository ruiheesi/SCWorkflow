#' PCA & Normalization template in single-cell-rna-seq-r4 NIDAP environment
#' from v96
#' 
#' @title PCA & Normalization  
#' @description This template normalizes each sample and plots Seurat PCAs and Elbow Plots for each sample to aid in selecting the number of PCs to carry forward in your analysis (see next template in workflow). PCA plots before and after regression will be shown for each variable chosen. The elbow plot will show reduction in variance (y-axis) as a function of PC count (x-axis). It is recommended to pick the number of PCs where variance is reduced by 98% (i.e. the first point where y coord = 0.02). This is Step 3 in the basic Single-Cell RNA-seq workflow.
#' @details
#' 
#'
#' @param Seurat_Object Please input a post-Filter QC Seurat Object to do your reduction on.
#' @param vars_to_regress Subtract (‘regress out’) this source of heterogeneity from the data. For example, to regress out mitochondrial effects, input "percent.mt." Options: percent.mt, nCount_RNA, S.Score, G2M.Score, CC.Difference
#' @param vars_to_plot Choose variables to visualize on PCA plot (maximum = 3)
#' @param npcs Select initial set of principal components
#' @param nfeatures Number of variable features
#' @param low_cut description.
#' @param high_cut description.
#' @param low_cut_disp description.
#' @param high_cut_disp description.
#' @param selection_method Method to choose top variable features. Options: vst, mean.var.plot, dispersion
#' @param doJackStraw Opt to visualize your data in a Jackstraw plot (sometimes more description than an elbow plot).
#' @param JackStraw_dims Recommended max 10
#' @param methods_PCA Methods available: Marchenko-Pastur: use eigenvalue null upper bound from URD, Elbow: Find threshold where percent change in variation between consecutive PCs is less than X% (set below)
#' @param var_threshold For Elbow method, set percent change threshold in variation between consecutive PCs
#' @param imageType Remember that svgs are much larger than pngs, so we recommend doing everything first in png, then rerunning to output specific svgs as needed. Options: png, svg
#' 
#' 
#' @import Seurat
#' @import gridExtra 
#' @import ggplot2
#' @import tidyverse
#' @import dplyr
#' @import svglite
#' @importFrom digest digest
#' @importFrom svglite svglite
#' 
#' 
#' @export
#' 
#' @return Seurat Objects and QC plots. This template normalizes each sample and plots Seurat PCAs and Elbow Plots for each sample to aid in selecting the number of PCs to carry forward in your analysis (see next template in workflow). PCA plots before and after regression will be shown for each variable chosen. The elbow plot will show reduction in variance (y-axis) as a function of PC count (x-axis). It is recommended to pick the number of PCs where variance is reduced by 98% (i.e. the first point where y coord = 0.02). This is Step 3 in the basic Single-Cell RNA-seq workflow.


PCA_and_Normalization <- function(Seurat_Object,
                           vars_to_regress = NULL,
                           vars_to_plot = c('percent.mt','Phase'),
                           npcs = 30,
                           nfeatures = 2000,
                           low_cut = 1,
                           high_cut = 8,
                           low_cut_disp = 1,
                           high_cut_disp = 100000,
                           selection_method = 'vst',
                           doJackStraw = FALSE,
                           JackStraw_dims = 5,
                           methods_PCA = c('Elbow','Marchenko-Pastur'),
                           var_threshold = 0.1,
                           imageType = 'png'
                           ){
    
  
 


  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  ################################
  # Cell Cycle Scoring and Find Variable Features

  CC_FVF_so <- function(so){
    # so <- CellCycleScoring(object = so, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
    # so$CC.Difference <- so$S.Score - so$G2M.Score
    so <- FindVariableFeatures(object = so, nfeatures = nfeatures, mean.cutoff = c(low_cut, high_cut), dispersion.cutoff = c(low_cut_disp, high_cut_disp), selection.method = selection_method)
    all.genes <- rownames(so)
    return(so)
  }
  
  ################################
  # Make PCA without regressing anything, and using only SCTransform().

  pca_noregress <- function(so) {
    so <- SCTransform(so,do.correct.umi = FALSE,return.only.var.genes = FALSE)
    so <- RunPCA(object = so, features = VariableFeatures(object = so), npcs = npcs)
    return(so)
  }
  
  ################################
  # Make PCA with SCTransform() (ScaleData deprecated) on all genes.
  pca <- function(so) {
    
    # Run SCTransform().        
    so <- SCTransform(so,do.correct.umi = TRUE, vars.to.regress = vars_to_regress, return.only.var.genes = FALSE)
    
    # Make PCA using last transform run, which will always be that from
    # SCTransform().
    so <- RunPCA(object = so, npcs = npcs)
    slot(so,"commands") <- list()
    return(so)
  }
  
  ################################
  ## Image Functions
  
  plotPCA <- function(so,m){
    p1 <- DimPlot(so, reduction = "pca")
    clusmat=data.frame(umap1=p1$data$PC_1,umap2=p1$data$PC_2, clusid=so@meta.data[[m]])
    
    sumpcsd = sum(so@reductions$pca@stdev)
    pcvar = (so@reductions$pca@stdev/sumpcsd)*100
    pcvar=formatC(pcvar,format = "g",digits=3)
    pcvar1 = pcvar[1] 
    pcvar2 = pcvar[2]
    pcvar 
    
    run.categ <- function(mat){
      colors=c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324","#800000","#aaffc3","#808000","#000075","#a9a9a9")
      
      g <- ggplot(mat, aes(x=umap1, y=umap2)) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        geom_point(aes(colour=clusid),size=0.5) +
        scale_color_manual(values=colors) +
        xlab(paste0("PC-1 ",pcvar[1],"%")) + ylab(paste0("PC-2 ",pcvar[2],"%")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),legend.text=element_text(size=rel(1.5))) +
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
        ggtitle(so@project.name)  
      return(g)
    }
    
    run.cont <- function(mat,midpt,maxpt){
      g <- ggplot(mat, aes(x=umap1, y=umap2)) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        geom_point(aes(colour=clusid),size=0.5) +
        scale_colour_gradient2(low = "blue",mid="lightgrey",high = "red",limits=c(0, maxpt),midpoint = midpt) +
        xlab(paste0("PC-1 ",pcvar[1],"%")) + ylab(paste0("PC-2 ",pcvar[2],"%")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),legend.text=element_text(size=rel(1.5))) +
        guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) +
        ggtitle(paste0(so@project.name,"_",m))  
      return(g)
      
      
    }
    if(class(clusmat$clusid) == "factor"){
      g = run.categ(clusmat)
      return(g)
    }
    else{
      clusmat %>% arrange(clusid) -> clusmat
      if(m=="percent.mt"){
        #mid = quantile(clusmat$clusid)[4]
        mid=5
        max = 10
        clusmat$clusid[clusmat$clusid > max] <- max
      }
      else{
        mid = quantile(clusmat$clusid)[3]
        max = quantile(clusmat$clusid,probs=0.95)
        clusmat$clusid[clusmat$clusid > max] <- max
      }
      
      g = run.cont(clusmat,mid,max)
      return(g)
    }
    
  }
  
  
  plotElbow <- function(so){
    
    if("Elbow" %in% methods_PCA){
      #Find Elbow:
      sumpcsd = sum(so@reductions$pca@stdev)
      pcvar = (so@reductions$pca@stdev/sumpcsd)*100
      cumu <- cumsum(pcvar)
      co1 <- which(cumu > 80 & pcvar < 5)[1]
      co2 <- sort(which((pcvar[1:length(pcvar) - 1] - pcvar[2:length(pcvar)]) > var_threshold), decreasing = T)[1] + 1
      pcs = min(co1,co2)
      lab = paste0("Elbow = ", pcs)
      xpos = pcs + 4
    }
    
    if("Marchenko-Pastur" %in% methods_PCA){
      #Using code from URD (https://rdrr.io/github/farrellja/URD/src/R/pca.R)
      pcaMarchenkoPastur <- function(M, N, pca.sdev, factor=1, do.print=T) {
        pca.eigenvalue <- (pca.sdev)^2
        marchenko.pastur.max <- (1+sqrt(M/N))^2
        pca.sig <- pca.eigenvalue > (marchenko.pastur.max * factor)
        if (do.print) {
          print(paste("Marchenko-Pastur eigenvalue null upper bound:", marchenko.pastur.max))
          if (factor != 1) {
            print(paste(length(which(pca.sig)), "PCs have eigenvalues larger than", factor, "times null upper bound."))
          } else {
            print(paste(length(which(pca.eigenvalue > marchenko.pastur.max)), "PCs have larger eigenvalues."))
          }}
        pca.sig.length = length(pca.sig[pca.sig==TRUE])
        return(pca.sig.length)
      }
      
      M = dim(so$RNA@data)[1]
      N = dim(so$RNA@data)[2]
      pca.sdev = so@reductions$pca@stdev
      pca.sig.num = pcaMarchenkoPastur(M=M,N=N,pca.sdev = pca.sdev)
      lab2 = paste0("MP = ", pca.sig.num)
      xpos2 = pca.sig.num+4
    }
    ep <- ElbowPlot(so,ndims=30) + theme_bw() + 
      ggtitle(paste0(so@project.name," Elbow Plot")) 
    
    if(exists("lab")){
      ep <- ep + 
        geom_vline(xintercept = pcs, color="red") +
        annotate("text",  x=xpos, y = 4, label = lab, color="red",size=4) 
    }
    if(exists("lab2")){
      ep <- ep + 
        geom_vline(xintercept = pca.sig.num, color="blue") +
        annotate("text",  x=xpos2, y = 6, label = lab2, color="blue",size=4)
    }
    return(ep)
    #print(ep)
    #return(NULL)
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  
  # # load data
  # object.class <- getClass(class(SO))
  # 
  # if(object.class@className == "RFoundryObject") {
  #   cat("1. Reading Seurat Object from dataset: RObjectdata.rds\n\n")
  #   
  #   SO = SO$value
  # } else {
  #   cat("1. Reading Seurat Object from dataset: SO.rds\n\n")
  #   
  #   fs <- SO$fileSystem()
  #   path <- fs$get_path("seurat_object.rds", 'r')
  #   SO <- readRDS(path)
  #   
  # }
  
  SO <- Seurat_Object

  
  #in case you want to redo this on a merged SO
  if (class(SO) =="Seurat") {
    x =list()
    x[[1]] <- SO
    SO <- x
  }
  
  ################################
  # Do transformation with and without regression using SCTransform().

  so_prep <- lapply(SO, CC_FVF_so) 
  so_orig <- lapply(so_prep, pca_noregress)
  so_list <- lapply(so_prep, pca) 
  
  
  ################################
  # Create Plots 
  
  imageCols = 2
  if (length(vars_to_plot) > 0) {
    imageCols <- imageCols + length(vars_to_plot)
  }   
  
  if(is.null(vars_to_plot)){
    vars_to_plot = "nCount_RNA"
  }
  len <- length(vars_to_plot)*2
  grobsList <- vector(mode = "list", length = len)
  
  
  
  ## PCA plots
  
  k=1

  for (i in 1:length(vars_to_plot)){ 
    grob <- lapply(so_orig, function(x) plotPCA(x,vars_to_plot[i]))
    grob=grid.arrange(grobs=grob,nrow=length(grob))
    grobsList[[k]] <- grob
    grob2 <- lapply(so_list, function(x) plotPCA(x,vars_to_plot[i]))
    grob2=grid.arrange(grobs=grob2,nrow=length(grob2))
    l=k+1
    grobsList[[l]] <- grob2
    k=k+2   
  }
  
  ## Elbow Plot
  

  grob3 <- lapply(so_list, function(x) plotElbow(x))
  grob3=grid.arrange(grobs=grob3,nrow=length(grob3))
  
  grobsList[[length(grobsList)+1]] <- grob3
  
  
  ## Jakstraw plot
  
  if (doJackStraw) {
    imageCols <- imageCols + 2
  }   
  if (doJackStraw) {
    grob4 <-lapply(so_list, function(x) JackStraw(x, reduction = "pca", dims = JackStraw_dims,
                                                  num.replicate = 100, prop.freq = 0.01, verbose = TRUE,
                                                  maxit = 1000))
    grob4 <- lapply(grob4, function(x) ScoreJackStraw(x, dims = 1:JackStraw_dims))
    grob4 <- lapply(grob4, function(x) JackStrawPlot(x, dims = 1:JackStraw_dims))
    grob4=grid.arrange(grobs=grob4,nrow=length(grob4))
    grobsList[[length(grobsList)+1]] <- grob4
  }
  
  grobs <- arrangeGrob(grobs=grobsList,ncol=length(grobsList),newpage=F)
  
  
  ################################
  # Create Figure Output 
  
  
  imageWidth = 1000*2*imageCols
  imageHeight = 1000*length(so_list)
  dpi = 300
  
  if (imageType == "png") {
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

  
  # plot(grobs)
  # grobs = arrangeGrob(grobs=grobs,nrow=1)
  # grid.arrange(grobs)
  grobs=grobs
  

  
  cat("\nPCA Object Checksum:\n")

  return(list(so=SO,plot=grobs))
  
  # output <- new.output()
  # output_fs <- output$fileSystem()
  # saveRDS(so_list, output_fs$get_path("seurat_object.rds", 'w'))
  # 
  # return(output_fs)
}
