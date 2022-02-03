library(png)
library(ggplot2)
library(grid)
library(reshape2)
library(spatstat.core)


setwd("~/Documents/main_files/AskExplain/generative_encoder/")


load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gcode___100.all.models.Rdata")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/BC_data_pixel.RData")

x_perturb <- BC_data_pixel

rm(BC_data_pixel)
gc()

epithelial_set <- "^WNT|^HSP|^FKBP|^AGRN|^AMBN|^AMLEX|^BMP|^BRCA|^COL|^EFE|^FBLN|^FBN|^IGF|^LTB|^MGP|^RSP|^SMOC|^TGFB|^ZP|CAN$|^VIM|^CD44|^ANXA1|^ACTA2|^ITGA8|^FN1|^VCAM1|^ITGB2|^CAV1|^LAM"
immune_set <- "^HLA|^CCR|^CCX|^CCL|^IL|^FCG|^GZM|^DUSP|^JUN|^NK|^CD1|^CD2|^CD3|^CD4|^CD5|^CD6|^CD7|^CD8|^CD9"

perturbation <- function(image_histology,
                         gcode_model_list,
                         gene_list,
                         perturb_gene_set, 
                         delta = "+",
                         factor_delta = 1){
  
  delta_image_histology <- 0
  for (gcode_model in gcode_model_list){
  gene_feature <- matrix(image_histology%*%gcode_model$main.parameters$beta[[2]]%*%t(gcode_model$main.parameters$beta[[1]]),nrow = 1)
  colnames(gene_feature) <- gene_list
  
  perturb_gene_ids <- gene_list[grep(gene_list, pattern = perturb_gene_set)]
  
  if (delta == "+"){
    gene_feature[,perturb_gene_ids] <- gene_feature[,perturb_gene_ids]+factor_delta*sd(gene_feature[,perturb_gene_ids])
  }
  if (delta == "-"){
    gene_feature[,perturb_gene_ids] <- gene_feature[,perturb_gene_ids]-factor_delta*sd(gene_feature[,perturb_gene_ids])
  }
  
  
  delta_image_histology <-  delta_image_histology + gene_feature%*%gcode_model$main.parameters$beta[[1]]%*%t((gcode_model$main.parameters$beta[[2]]))
  }
  
  return(delta_image_histology)
}





convert_to_RGB <- function(x){
  e <- ecdf(x)
  j <- e(x)
  x <- array(j,dim(x))
  return(x)
}






for (image_id in c(1:250)){
  
  gg_plots <- list()
  
  
  for (delta_change in list(c("( - ) epithelial","( - ) immune"),"( - ) epithelial",c("( - ) epithelial","( + ) immune"),"( - ) immune","Original tissue","( + ) immune",c("( + ) epithelial","( - ) immune"),"( + ) epithelial",c("( + ) epithelial","( + ) immune"))){
    
    x_perturb.copy <- x_perturb[image_id,,drop=F]
    
    if ("( + ) epithelial" %in% delta_change){
      
      x_perturb.copy <- perturbation(image_histology = x_perturb.copy, 
                                     gcode_model = gcode.all.models$gcode.non_tumour, 
                                     gene_list = gene_consensus$hgnc_symbol,
                                     perturb_gene_set = epithelial_set, 
                                     delta = "+", 
                                     factor_delta = 10)

    }
    
    
    if ("( - ) epithelial" %in% delta_change){
      
      x_perturb.copy <- perturbation(image_histology = x_perturb.copy, 
                                     gcode_model = gcode.all.models$gcode.non_tumour, 
                                     gene_list = gene_consensus$hgnc_symbol,
                                     perturb_gene_set = epithelial_set, 
                                     delta = "-", 
                                     factor_delta = 10)
      
    }
    
    
    if ("( + ) immune" %in% delta_change){
      
      x_perturb.copy <- perturbation(image_histology = x_perturb.copy, 
                                     gcode_model = gcode.all.models$gcode.non_tumour,
                                     gene_list = gene_consensus$hgnc_symbol,
                                     perturb_gene_set = immune_set, 
                                     delta = "+", 
                                     factor_delta = 10)
      
    }
    
    
    if ("( - ) immune" %in% delta_change){
      
      x_perturb.copy <- perturbation(image_histology = x_perturb.copy, 
                                     gcode_model = gcode.all.models$gcode.non_tumour, 
                                     gene_list = gene_consensus$hgnc_symbol,
                                     perturb_gene_set = immune_set, 
                                     delta = "-", 
                                     factor_delta = 10)
      
    }
    
    
    g <- (convert_to_RGB(aperm(array((x_perturb.copy),dim=c(64,64,3)),c(1,2,3))))
    g <- rasterGrob(g, interpolate=TRUE)
    g_plots <- ggplot2::qplot(1:10, 1:10, geom="blank") +
      annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
      ggtitle(paste0(delta_change,collapse="  ,  "))
    g_plots
    gg_plots <- c(gg_plots,list(g_plots))
    
  }
  
  
  
  
  gg_plots <- lapply(gg_plots,function(X){
    X <- X + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                   axis.text.y=element_blank(),axis.ticks=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),legend.position="none",
                   panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),plot.background=element_blank())
    X <- X + theme(plot.title = element_text(size = 16, face = "bold", color="black", hjust=0.5))
  })
  
  
  
  lm <- rbind(c(1,2,3),
              c(4,5,6),
              c(7,8,9))
  
  
  library(grid)
  library(gridExtra)
  final_plots <- arrangeGrob(
    grobs = gg_plots,
    layout_matrix = lm
  )
  
  plot(final_plots)
  ggsave(filename = paste("./figures/explore/jpeg_deep_immune_epithelial/",image_id,"___epithelial_vertical___immune_horizontal___transition.png",sep=""),plot = final_plots, width = 10, height = 10)
  
  
  
}

