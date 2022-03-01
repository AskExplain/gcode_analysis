setwd("~/Documents/main_files/AskExplain/generative_encoder/")
library(lsa)

load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/spinal/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/spinal/gcode___spinal.5.all.models.Rdata")

hires_image <- list.files("./data/external/spinal/images/HE/",recursive = T,pattern = "*.jpg")
coordinate_list <- list.files("./data/external/spinal/annotations/",recursive = T,pattern = "*.tsv")
gex_raw <- list.files("./data/external/spinal/count_matrices/",recursive = T,pattern = "*.txt")

image_IDS <- do.call('c',lapply(strsplit(coordinate_list,split = "\\."),function(X){X[1]}))

main_path <- c("~/Documents/main_files/AskExplain/generative_encoder/data/external/spinal/")

extract_pixels <- function(image,coords){
  
  count <- 1
  coord_list <- list()
  for (i in c(0,1)){
    for (j in c(0,1)){
      copy_coords <- coords
      copy_coords[,1] <- coords[,1] - 32/2 + i * 32/2 
      copy_coords[,2] <- coords[,2] - 32/2 + j * 32/2
      
      coord_list[[count]] <- do.call('rbind',pbmcapply::pbmclapply(c(1:dim(copy_coords)[1]),function(i){
        
        cropped_image <- magick::image_crop(image = image, geometry = paste("64x64+",copy_coords[i,1],"+",copy_coords[i,2],sep=""))
        
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


set.seed(1)
train_IDS <- sample(c(1:length(image_IDS)),length(image_IDS)*0.35)
test_IDS <- sample(c(1:length(image_IDS))[-train_IDS],100)



spinal_validation <- list(validation=list(NULL),gex=list(NULL),pixel=list(NULL))


for (image_id in c(1:length(test_IDS))){
  
  spinal_validation$gex[[image_id]] <- t(read.delim2(paste(main_path,"/count_matrices/",gex_raw[test_IDS[image_id]],sep=""),header=T,sep="\t",row.names=NULL))
  colnames(spinal_validation$gex[[image_id]]) <- spinal_validation$gex[[image_id]][1,]
  row_ids <- row.names(spinal_validation$gex[[image_id]])[-1]
  spinal_validation$gex[[image_id]] <- apply(spinal_validation$gex[[image_id]][-1,],2,as.numeric)
  gene_ids <- colnames(spinal_validation$gex[[image_id]])
  
  observed_gex <- spinal_validation$gex[[image_id]][,gene_consensus]
  spinal_validation$gex[[image_id]] <- NULL
  
  coordinate_file <- apply(do.call('rbind',strsplit(gsub(row_ids,pattern = "X",replacement = ""),split = "_")),2,as.numeric)
  
  histo <- magick::image_read(paste(main_path,"/images/HE/",hires_image[test_IDS[image_id]],sep=""))
  main_coords <- apply(as.matrix(apply(coordinate_file,2,as.numeric))*194/6400,2,function(X){pmin(X,1)})*as.numeric(magick::image_info(histo)[2:3])
  main_image <- extract_pixels(image = histo, coords = main_coords)  
  
  if (!any(do.call('c',lapply(c(1:4),function(X){is.character(main_image[[X]])})))){

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
    
  } else {
    to_return <- NULL
  }
  
  spinal_validation$pixel[[image_id]] <- NULL
  spinal_validation$gex[[image_id]] <- NULL
  
  spinal_validation$validation[[image_id]] <- to_return
  
  
}


save(spinal_validation, file = "./data/workflow/spinal/spinal_validation.RData")







load("./data/workflow/spinal/spinal_validation.RData")
source("./github/accuracy_gcode_spinal_SPATIAL/3_evaluate/plot_similarity_predict_vs_observe.R")


gene_similarity <- do.call('rbind',lapply(spinal_validation$validation,function(list_quadrants){
  
  do.call('rbind',lapply(list_quadrants,function(X){
    X$gene
  }))
  
}))


sample_similarity <- do.call('rbind',lapply(spinal_validation$validation,function(list_quadrants){
  
  do.call('rbind',lapply(list_quadrants,function(X){
    X$sample
  }))
  
}))


plot_similarity_predict_vs_observe(sample_similarity = sample_similarity, feature_similarity = gene_similarity, title_name = "spatial_perturbation", tissue_name = "spinal")





