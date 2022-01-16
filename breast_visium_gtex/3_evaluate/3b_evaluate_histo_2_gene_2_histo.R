library(png)
library(ggplot2)
library(grid)
library(reshape2)
library(spatstat.core)



convert_to_RGB <- function(x){
  e <- ecdf(x)
  j <- e(x)
  x <- array(j,dim(x))
  return(x)
}


setwd("~/Documents/main_files/AskExplain/generative_encoder/")


load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gcode.all.models.Rdata")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/BC_data_pixel.RData")

x <- BC_data_pixel

rm(BC_data_pixel)
gc()

x_perturb <- x%*%gcode.all.models$gcode.non_tumour$main.parameters$beta[[2]]%*%t(gcode.all.models$gcode.non_tumour$main.parameters$beta[[1]])
colnames(x_perturb) <- gene_consensus$hgnc_symbol


for (image_id in c(1:250)){
  
  
  
  gg_plots <- list()
  for (delta_change in list(c("( - ) epithelial","( - ) immune"),"( - ) epithelial",c("( - ) epithelial","( + ) immune"),"( - ) immune","Original tissue","( + ) immune",c("( + ) epithelial","( - ) immune"),"( + ) epithelial",c("( + ) epithelial","( + ) immune"))){
    
    x_perturb.copy <- x_perturb
    
    if ("( + ) epithelial" %in% delta_change){
      
      epithelial_gene_ids <- gene_consensus$hgnc_symbol[grep(gene_consensus$hgnc_symbol, pattern = "^WNT|^HSP|^FKBP|^AGRN|^AMBN|^AMLEX|^BMP|^BRCA|^COL|^EFE|^FBLN|^FBN|^IGF|^LTB|^MGP|^RSP|^SMOC|^TGFB|^ZP|CAN$|^VIM|^CD44|^ANXA1|^ACTA2|^ITGA8|^FN1|^VCAM1|^ITGB2|^CAV1|^LAM")]
      x_perturb.copy[,epithelial_gene_ids] <- x_perturb.copy[,epithelial_gene_ids]+qnorm(0.125)*sd(x_perturb.copy[,epithelial_gene_ids])
      
    }
    
    
    if ("( - ) epithelial" %in% delta_change){
      
      epithelial_gene_ids <- gene_consensus$hgnc_symbol[grep(gene_consensus$hgnc_symbol, pattern = "^WNT|^HSP|^FKBP|^AGRN|^AMBN|^AMLEX|^BMP|^BRCA|^COL|^EFE|^FBLN|^FBN|^IGF|^LTB|^MGP|^RSP|^SMOC|^TGFB|^ZP|CAN$|^VIM|^CD44|^ANXA1|^ACTA2|^ITGA8|^FN1|^VCAM1|^ITGB2|^CAV1|^LAM")]
      x_perturb.copy[,epithelial_gene_ids] <- x_perturb.copy[,epithelial_gene_ids]-qnorm(0.125)*sd(x_perturb.copy[,epithelial_gene_ids])      
    
    }
    
    
    if ("( + ) immune" %in% delta_change){
      
      immune_gene_ids <- gene_consensus$hgnc_symbol[grep(gene_consensus$hgnc_symbol, pattern = "^HLA|^CCR|^CCX|^CCL|^IL|^FCG|^GZM|^DUSP|^JUN|^NK|^CD1|^CD2|^CD3|^CD4|^CD5|^CD6|^CD7|^CD8|^CD9")]
      x_perturb.copy[,immune_gene_ids] <- x_perturb.copy[,immune_gene_ids]+qnorm(0.125)*sd(x_perturb.copy[,immune_gene_ids])      
    
    }
    
    
    if ("( - ) immune" %in% delta_change){
      
      immune_gene_ids <- gene_consensus$hgnc_symbol[grep(gene_consensus$hgnc_symbol, pattern = "^HLA|^CCR|^CCX|^CCL|^IL|^FCG|^GZM|^DUSP|^JUN|^NK|^CD1|^CD2|^CD3|^CD4|^CD5|^CD6|^CD7|^CD8|^CD9")]
      x_perturb.copy[,immune_gene_ids] <- x_perturb.copy[,immune_gene_ids]-qnorm(0.125)*sd(x_perturb.copy[,immune_gene_ids])      
      
    }
    
    
    
    g <- (convert_to_RGB(aperm(array((x_perturb.copy[image_id,,drop=F]%*%(gcode.all.models$gcode.non_tumour$main.parameters$beta[[1]])%*%t(gcode.all.models$gcode.non_tumour$main.parameters$beta[[2]])),dim=c(100,100,3)),c(1,2,3))))
    g <- rasterGrob(g, interpolate=TRUE)
    g_plots <- ggplot2::qplot(1:10, 1:10, geom="blank") +
      annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
      ggtitle(paste0(delta_change,collapse="  ,  "))
    
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
  
  ggsave(filename = paste("./figures/explore/epithelial_immune_transition/",image_id,"___epithelial_vertical___immune_horizontal___transition.png",sep=""),plot = final_plots, width = 10, height = 10)
  
  
  
}




write.csv(paste0(sort(gene_consensus$hgnc_symbol[grep(gene_consensus$hgnc_symbol, pattern = "^WNT|^HSP|^FKBP|^AGRN|^AMBN|^AMLEX|^BMP|^BRCA|^COL|^EFE|^FBLN|^FBN|^IGF|^LTB|^MGP|^RSP|^SMOC|^TGFB|^ZP|CAN$|^VIM|^CD44|^ANXA1|^ACTA2|^ITGA8|^FN1|^VCAM1|^ITGB2|^CAV1|^LAM")]),collapse=", "),file = "./data/workflow/epithelial_genes.csv",quote=F,row.names = F)
write.csv(paste0(sort(gene_consensus$hgnc_symbol[grep(gene_consensus$hgnc_symbol, pattern = "^HLA|^CCR|^CCX|^CCL|^IL|^FCG|^GZM|^DUSP|^JUN|^NK|^CD1|^CD2|^CD3|^CD4|^CD5|^CD6|^CD7|^CD8|^CD9")]),collapse=", "),file = "./data/workflow/immune_genes.csv",quote=F,row.names = F)

