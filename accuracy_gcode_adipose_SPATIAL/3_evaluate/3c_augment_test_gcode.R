
setwd("~/Documents/main_files/AskExplain/generative_encoder/")
library(lsa)


load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/gcode___adipose.5.all.models.Rdata")


hires_image <- list.files("./data/external/adipocyte//",recursive = T,pattern = "*hires_image.png")
scaling_json <- list.files("./data/external/adipocyte//",recursive = T,pattern = "*json.json")
coordinate_list <- list.files("./data/external/adipocyte//",recursive = T,pattern = "*list.csv")
gex_raw <- list.files("./data/external/adipocyte//",recursive = T,pattern = "*.mtx.gz")
gex_raw <- gex_raw[grep("raw",gex_raw)]
gex_features <- list.files("./data/external/adipocyte//",recursive = T,pattern = "*features.tsv.gz")
gex_features <- gex_features[grep("raw",gex_features)]

image_IDS <- do.call('c',lapply(strsplit(coordinate_list,split = "/"),function(X){X[1]}))

main_path <- c("~/Documents/main_files/AskExplain/generative_encoder/data/external/adipocyte/")

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


set.seed(1)
train_IDS <- sample(c(1:length(image_IDS)),length(image_IDS)*0.6)
test_IDS <- c(1:length(image_IDS))[-train_IDS]

adipose_validation <- list(gex = list(NULL), pixel = list(NULL), validation = list(NULL))

for (image_id in c(1:length(test_IDS))){
  
  coordinate_file <- read.csv(paste(main_path,coordinate_list[test_IDS[image_id]],sep=""),header=F)
  scaling_file <- read.delim(paste(main_path,scaling_json[test_IDS[image_id]],sep=""),header=F)
  hires_scaling_factor <- as.numeric(gsub(" tissue_hires_scalef: ","",strsplit(scaling_file[1,1],split = ",")[[1]][2]))
  coordinate_file[,c(5,6)] <- coordinate_file[,c(5,6)]*hires_scaling_factor
  
  adipose_validation$gex[[image_id]] <- (Matrix::readMM(paste(main_path,gex_raw[test_IDS[image_id]],sep="")))
  gene_ids <- read.table(paste(main_path,gex_features[test_IDS[image_id]],sep=""))
  row.names(adipose_validation$gex[[image_id]]) <- gene_ids$V2
  
  IDS_remove <- which(Matrix::colSums(adipose_validation$gex[[image_id]])<2000)
  
  observed_gex <- t(as.matrix(adipose_validation$gex[[image_id]][,-IDS_remove]))
  histo <- magick::image_read(paste(main_path,hires_image[test_IDS[image_id]],sep=""))
  
  main_image <- extract_pixels(image = histo, coords = coordinate_file[-IDS_remove,c(5,6)])
  
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
  
  adipose_validation$validation[[image_id]] <- to_return
  
}




save(adipose_validation, file = "./data/workflow/adipose/adipose_validation.RData")









load("./data/workflow/adipose/adipose_validation.RData")
source("./github/accuracy_gcode_adipose_SPATIAL/3_evaluate/plot_similarity_predict_vs_observe.R")


gene_similarity <- do.call('rbind',lapply(adipose_validation$validation,function(list_quadrants){
  
  do.call('rbind',lapply(list_quadrants,function(X){
    X$gene
  }))
  
}))


sample_similarity <- do.call('rbind',lapply(adipose_validation$validation,function(list_quadrants){
  
  do.call('rbind',lapply(list_quadrants,function(X){
    X$sample
  }))
  
}))


plot_similarity_predict_vs_observe(sample_similarity = sample_similarity, feature_similarity = gene_similarity, title_name = "spatial_perturbation", tissue_name = "adipose")





