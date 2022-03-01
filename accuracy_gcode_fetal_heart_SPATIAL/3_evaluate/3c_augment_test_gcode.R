setwd("~/Documents/main_files/AskExplain/generative_encoder/")

library(lsa)

load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/fetal_heart/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/fetal_heart/gcode___fetal_heart.5.all.models.Rdata")

hires_image <- list.files("./data/external/fetal_heart/images/",recursive = T,pattern = "*HE_small.jpg")
coordinate_list <- list.files("./data/external/fetal_heart/images/",recursive = T,pattern = "*.tsv")
gex_raw <- c("./data/external/fetal_heart/matrix/Filtered/filtered_ST_matrix_and_meta_data/filtered_matrix.tsv")
meta_raw <- c("./data/external/fetal_heart/matrix/Filtered/filtered_ST_matrix_and_meta_data/meta_data.tsv")

image_IDS <- do.call('c',lapply(strsplit(coordinate_list,split = "/"),function(X){X[2]}))

main_path <- c("~/Documents/main_files/AskExplain/generative_encoder/data/external/fetal_heart/")

extract_pixels <- function(image,coords){
  
  count <- 1
  coord_list <- list()
  for (i in c(0,1)){
    for (j in c(0,1)){
      copy_coords <- coords
      copy_coords[,1] <- coords[,1] - 50 + i * 50
      copy_coords[,2] <- coords[,2] - 50 + j * 50
      
      coord_list[[count]] <- do.call('rbind',pbmcapply::pbmclapply(c(1:dim(copy_coords)[1]),function(i){
        
        cropped_image <- magick::image_crop(image = image, geometry = paste("100x100+",copy_coords[i,1],"+",copy_coords[i,2],sep=""))
        
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
train_IDS <- sample(c(1:length(image_IDS)),length(image_IDS)*0.6)
test_IDS <- c(1:length(image_IDS))[-train_IDS]

meta_table <- read.delim2(meta_raw,header=T,sep="\t",row.names=NULL)
gex_table <- read.delim2(gex_raw,header=T,sep="\t",row.names=NULL)
gene_ids <- gex_table[,1]
cell_ids <- colnames(gex_table)[-1]
gex_table <- gex_table[,-1]

fetal_heart_validation <- list(validation = list(NULL))

for (image_id in c(1:length(test_IDS))){
  
  coordinates <- gsub(meta_table[,1],pattern = "X",replacement = "")
  coordinate_file <- read.delim2(paste(main_path,"/images/",coordinate_list[image_id],sep=""),row.names=NULL)
  coordinate_file$coord <- paste(image_id,coordinate_file$x,coordinate_file$y,sep="x")
  
  coordinate_file <- coordinate_file[match(coordinates[grep(pattern = paste("^",image_id,"x",sep=""),coordinates)],coordinate_file$coord),]
  
  histo <- magick::image_read(paste(main_path,"/images/",hires_image[test_IDS[image_id]],sep=""))
  main_image <- extract_pixels(image = histo, coords = apply(coordinate_file[,c(5,6)],2,as.numeric))  
  observed_gex <- t(gex_table[,paste("X",coordinates[grep(pattern = paste("^",image_id,"x",sep=""),coordinates)],sep="")])
 
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
    
  }
  
  fetal_heart_validation$validation[[image_id]] <- to_return
  
  
}



save(fetal_heart_validation, file = "./data/workflow/fetal_heart/fetal_heart_validation.RData")



load("./data/workflow/fetal_heart/fetal_heart_validation.RData")
source("./github/accuracy_gcode_fetal_heart_SPATIAL/3_evaluate/plot_similarity_predict_vs_observe.R")


gene_similarity <- do.call('rbind',lapply(fetal_heart_validation$validation,function(list_quadrants){
  
  do.call('rbind',lapply(list_quadrants,function(X){
    X$gene
  }))
  
}))


sample_similarity <- do.call('rbind',lapply(fetal_heart_validation$validation,function(list_quadrants){
  
  do.call('rbind',lapply(list_quadrants,function(X){
    X$sample
  }))
  
}))


plot_similarity_predict_vs_observe(sample_similarity = sample_similarity, feature_similarity = gene_similarity, title_name = "spatial_perturbation", tissue_name = "fetal_heart")





