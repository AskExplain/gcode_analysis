
######## 	 All files can be found below:

### Spatial Transcriptomics dataset
# Data downloaded at 
# https://data.mendeley.com/datasets/29ntw7sh4r/5 and put in a file at ./data/external/breast/
###

###
# The R files for Generative Encoding can be found at
# https://drive.google.com/drive/folders/1a6txh_ZBcpPnLFj9mTknfwZ3A0dNTIlS
###

#########

library(lsa)

setwd("~/Documents/main_files/AskExplain/generative_encoder/")
ids <- list.files("./data/external/breast_visium/")
breast_cancer_ids <- substr(do.call('c',lapply(strsplit(ids[grep("Coords",ids)],split = "_"),function(X){paste(X[1],X[2],sep="_")})),3,10)


load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/breast/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/breast/gcode___breast.5.all.models.Rdata")


extract_pixels <- function(image,coords){
  
  count <- 1
  coord_list <- list()
  for (i in c(0,1)){
    for (j in c(0,1)){
      copy_coords <- coords
      copy_coords[,1] <- coords[,1] - 64 + i * 64 
      copy_coords[,2] <- coords[,2] - 64 + j * 64
      
      coord_list[[count]] <- do.call('rbind',pbmcapply::pbmclapply(c(1:dim(copy_coords)[1]),function(i){
        
        cropped_image <- magick::image_crop(image = image, geometry = paste("128x128+",copy_coords[i,1],"+",copy_coords[i,2],sep=""))
        
        coord_pixels <- c(as.numeric(magick::image_data(cropped_image, 'rgb')))
        
      },mc.cores = 4))
      
      count <- count + 1
    }
  }
  
  
  return(coord_list) 
}




fa0 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$beta[[1]]
fb0 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$beta[[2]]
fi0 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$intercept[[2]]
fi1 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$intercept[[1]]
fc <- gcode.all.models$gcode.non_tumour[[1]]$main.code$code[[1]]


breast_validation <- list(gex = list(NULL), pixel = list(NULL), validation = list(NULL))

for (id in breast_cancer_ids){
  
  print(id)
  
  BC_covar <- list()
  
  all_files <- ids[grep(id,ids)]
  
  coords <- read.table(paste("./data/external/breast_visium//",all_files[1],sep=""))
  gex <- read.table(paste("./data/external/breast_visium/",all_files[2],sep=""))
  histo <- magick::image_read(paste("./data/external/breast_visium/",all_files[3],sep=""))
  spots <- read.csv(paste("./data/external/breast_visium/",all_files[4],sep=""))
  
  spots_internal <- cbind(spots[spots[,1]%in%row.names(gex),],id)
  BC_covar$non_tumour <- spots_internal[coords$tumor=="non",]
  
  if (dim(BC_covar$non_tumour)[1]>0){
    observed_gex <- gex[row.names(gex)%in%BC_covar$non_tumour[,1],gene_consensus$ensembl_gene_id]
    main_image <- extract_pixels(image = histo, coords = BC_covar$non_tumour[,c(2,3)])

    to_return <- lapply(c(1:4),function(X){
      print(X)
      
      predicted_gex.alpha <- as.matrix((main_image[[X]] ) - c(fi0))%*%fb0%*%t(fc)%*%MASS::ginv(fc%*%t(fc))
      predicted_gex <- (predicted_gex.alpha%*%fc%*%t(fa0)+c(fi1))
      
      gene_main_metrics <- do.call('rbind',pbmcapply::pbmclapply(c(1:dim(observed_gex)[2]),function(id){
        
        main_cor <- c(cosine(as.numeric(observed_gex[,id]),as.numeric(predicted_gex[,id])))
        main_mae <- c(mean(abs(as.numeric(observed_gex[,id])-as.numeric(predicted_gex[,id]))))
        
        list(main_cor = main_cor,main_mae = main_mae)
      },mc.cores = 8))
      
      samples_main_metrics <- do.call('rbind',pbmcapply::pbmclapply(c(1:dim(observed_gex)[1]),function(id){
        
        main_cor <- c(cosine(as.numeric(observed_gex[id,]),as.numeric(predicted_gex[id,])))
        main_mae <- c(mean(abs(as.numeric(observed_gex[id,])-as.numeric(predicted_gex[id,]))))
        
        list(main_cor = main_cor,main_mae = main_mae)
      },mc.cores = 8))
      
      to_return <- list(gene=gene_main_metrics,sample=samples_main_metrics)
      
      return(to_return)
      
    })
    
    
    
    breast_validation$validation[[id]] <- to_return
  }

}



save(breast_validation, file = "./data/workflow/breast/breast_validation.RData")




load("./data/workflow/breast/breast_validation.RData")
source("./github/accuracy_gcode_breast_SPATIAL/3_evaluate/plot_similarity_predict_vs_observe.R")


gene_similarity <- do.call('rbind',lapply(breast_validation$validation,function(list_quadrants){
  
  do.call('rbind',lapply(list_quadrants,function(X){
    X$gene
  }))
  
}))


sample_similarity <- do.call('rbind',lapply(breast_validation$validation,function(list_quadrants){
  
  do.call('rbind',lapply(list_quadrants,function(X){
    X$sample
  }))
  
}))


plot_similarity_predict_vs_observe(sample_similarity = sample_similarity, feature_similarity = gene_similarity, title_name = "spatial_perturbation", tissue_name = "breast")





