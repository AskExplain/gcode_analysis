perturbation <- function(image_histology,
                         model_ensemble,
                         gene_list,
                         perturb_gene_set, 
                         delta = 1,
                         factor_delta = 0.01){
  
  ids <- grep(gene_list, pattern = perturb_gene_set)
  perturb_gene_ids <- gene_list[ids]
  
  gene_decoder <- gcode_model$main.parameters$beta[[1]]%*%t(gcode_model$main.parameters$beta.code[[1]])
  pixel_encoder <- (gcode_model$main.parameters$beta.code[[2]])%*%tgcode_model$main.parameters$beta[[2]]
  
  gene_feature <- matrix(image_histology%*%pixel_encoder%*%t(gene_decoder),nrow = 1)
  colnames(gene_feature) <- gene_list
  
  gene_feature[,perturb_gene_ids] <- gene_feature[,perturb_gene_ids]+delta*1.96*sd(gene_feature[,perturb_gene_ids])
  
  delta_image_histology <- gene_feature%*%((gene_decoder)%*%t(pixel_encoder))
  
  return(delta_image_histology)
}



convert_to_RGB <- function(x){
  e <- ecdf(x)
  j <- e(x)
  x <- array(j,dim(x))
  return(x)
}


plot_clean_ggplot2 <- function(image_list){
  
  character_main = c("(-) epithelial","Original tissue", "(+) epithelial")
  
  gg_plots <- lapply(c(1:length(image_list)),function(X){
    
    image_main = image_list[[X]]
    
    ggplot2::qplot(1:10, 1:10, geom="blank") +
      annotation_custom(
        
        rasterGrob((convert_to_RGB(aperm(array((image_main),dim=c(100,100,3)),c(1,2,3)))), interpolate=TRUE)
        
        , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
      ggtitle(paste0(character_main[X],collapse="  ,  "))
    
  })
  
  
  gg_plots <- lapply(gg_plots,function(X){
    X <- X + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                   axis.text.y=element_blank(),axis.ticks=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),legend.position="none",
                   panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),plot.background=element_blank())
    X <- X + theme(plot.title = element_text(size = 16, face = "bold", color="black", hjust=0.5))
  })
  
  
  
  lm <- rbind(c(1,2,3))
  
  final_plots <- arrangeGrob(
    grobs = gg_plots,
    layout_matrix = lm
  )
  
  return(final_plots)
}

