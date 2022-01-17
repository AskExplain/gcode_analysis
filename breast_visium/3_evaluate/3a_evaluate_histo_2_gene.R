load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gcode.all.models.Rdata")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gene_consensus.RData")


setwd("~/Documents/main_files/AskExplain/generative_encoder/")
ids <- list.files("./data/external/breast_visium/")
breast_cancer_ids <- substr(do.call('c',lapply(strsplit(ids[grep("Coords",ids)],split = "_"),function(X){paste(X[1],X[2],sep="_")})),3,10)



convert_to_HEAT <- function(x){
  x[x<0] <- x[x<0]*(-1) / min(x)
  x[x>0] <- x[x>0] / max(x)
  return(x)
}


convert_to_RGB <- function(x){
  e <- ecdf(x)
  j <- e(x)
  x <- array(j,dim(x))
  return(x)
}


extract_pixels <- function(image,coords){
  
  coord_pixels <- c()
  
  for (i in 1:dim(coords)[1]){
    
    if (i%%100 == 0){
      print(i)
    }
    
    X <- coords[i,1] 
    Y <- coords[i,2] 
    
    cropped_image <- magick::image_crop(image = image, geometry = paste("100x100+",X,"+",Y,sep=""))
    coord_pixels <- rbind(coord_pixels,c(as.numeric(magick::image_data(cropped_image, 'rgb'))))
  }
  
  return(coord_pixels) 
}


library(png)
library(ggplot2)
library(grid)
library(reshape2)
library(spatstat.core)


for (id in breast_cancer_ids){
  
  print(id)
  
  all_files <- ids[grep(id,ids)]
  
  
  histo <- magick::image_read(paste("./data/external/breast_visium/",all_files[3],sep=""))
  histo_gg2 <- magick::image_ggplot(histo)
  
  step <- 50
  histo_image <- magick::image_data(histo, 'rgb')
  dimensions_full <- dim(histo_image)
  dimensions <- dimensions_full[-1] - (100+1)
  histo_image <- magick::image_crop(image = histo, geometry = paste(floor(dimensions[1]/step)*step,"x",floor(dimensions[2]/step)*step,"+0+0",sep=""))
  histo_numeric <- magick::image_data(histo_image, 'rgb')
  dimensions_final <- dim(histo_numeric)
  
  rm(histo)
  rm(histo_image)
  gc()
  
  coords <- cbind(rep(seq((0),(dimensions_final[2]-2*step),step),each=length(seq((0),(dimensions_final[3]-2*step),step))),rep(seq((0),(dimensions_final[3]-2*step),step),times=length(seq((0),(dimensions_final[2]-2*step),step))))
  
  print(dimensions_final)
  print(apply(coords,2,max))
  print(dim(coords))
  
  histo_pixel_matrix <- do.call('rbind',pbmcapply::pbmclapply(c(1:dim(coords)[1]),function(X){
    c(as.numeric(histo_numeric[,(coords[X,1]+1):(coords[X,1]+100),(coords[X,2]+1):(coords[X,2]+100)]))%*%gcode.all.models$gcode.non_tumour$main.parameters$beta[[2]]
  },mc.cores = 8))
  
  histo_main_matrix <- histo_pixel_matrix%*%t(gcode.all.models$gcode.non_tumour$main.parameters$beta[[1]])

  colnames(histo_main_matrix) <- gene_consensus$hgnc_symbol
  
  # gene_consensus$hgnc_symbol[grep(x = gene_consensus$hgnc_symbol,pattern = "TNFS")]
  
  g_list <- list()
  for (gene_of_interest in c("DNASE1","RNASE1","BRCA1",
                             "PRF1","IL2RB","NKG7",
                             "BAX","BCL2","BCL3",
                             "TLR2","TLR4","TLR5",
                             "TNFSF10","TNFSF12","TNFSF13B",
                             "GNAS","FASN","DDX5",
                             "BGN","ACTG1","AEBP1")){
    print(gene_of_interest)
    g <- melt(as.matrix(as.im((apply(t(matrix(((histo_main_matrix[,gene_of_interest])),nrow=dimensions_final[2]/step-1,ncol=dimensions_final[3]/step-1,byrow = T)),2,rev)))))
    g$value <- convert_to_HEAT(scale(g$value))
    heatmap_p <- ggplot2::ggplot(data=g, aes(x=Var2, y=Var1, fill=value))+ 
      geom_tile(color="white") +
      scale_fill_gradient2(low = "white", high = "red",midpoint=0,limit = c(-1,1)) +
      ggtitle(gene_of_interest)
    plot(heatmap_p)
    g_list <- c(g_list,list(heatmap_p))
    
  }
  
  g_list <- c(list(histo_gg2),g_list)
  
  
  g_list <- lapply(g_list,function(X){
    X <- X + theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank())
    
    X <- X + theme(plot.title = element_text(size = 40, face = "bold"))
  })
  
  
  
  lm <- rbind(cbind(rbind(c(1,1),c(1,1)),rbind(c(2,3,4),c(5,6,7))),rbind(c(8,11,14,17,20),c(9,12,15,18,21),c(10,13,16,19,22)))
  
  library(grid)
  library(gridExtra)
  g_plots <- arrangeGrob(
    grobs = g_list,
    layout_matrix = lm
  )
  

  ggsave(filename = paste("./figures/explore/heatmap_he_et_al/",id,"___he_et_al_ST_breast_tissue.png",sep=""),plot = g_plots, width = 30,height = 30)
  
}

