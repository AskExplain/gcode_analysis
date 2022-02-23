setwd("~/Documents/main_files/AskExplain/generative_encoder/")

hires_image <- list.files("./data/external/fetal_heart/images/",recursive = T,pattern = "*HE_small.jpg")
coordinate_list <- list.files("./data/external/fetal_heart/images/",recursive = T,pattern = "*.tsv")
gex_raw <- c("./data/external/fetal_heart/matrix/Filtered/filtered_ST_matrix_and_meta_data/filtered_matrix.tsv")
meta_raw <- c("./data/external/fetal_heart/matrix/Filtered/filtered_ST_matrix_and_meta_data/meta_data.tsv")

image_IDS <- do.call('c',lapply(strsplit(coordinate_list,split = "/"),function(X){X[2]}))

main_path <- c("~/Documents/main_files/AskExplain/generative_encoder/data/external/fetal_heart/")

extract_pixels <- function(image,coords){
  
  coord_pixels <- do.call('rbind',pbmcapply::pbmclapply(c(1:dim(coords)[1]),function(i){
    X <- coords[i,1] - 50
    Y <- coords[i,2] - 50
    
    cropped_image <- magick::image_crop(image = image, geometry = paste("100x100+",X,"+",Y,sep=""))
    coord_pixels <- c(as.numeric(magick::image_data(cropped_image, 'rgb')))
    
  },mc.cores = 8))
  
  return(coord_pixels) 
}

set.seed(1)
train_IDS <- sample(c(1:length(image_IDS)),length(image_IDS)*0.6)
test_IDS <- c(1:length(image_IDS))[-train_IDS]

meta_table <- read.delim2(meta_raw,header=T,sep="\t",row.names=NULL)
gex_table <- read.delim2(gex_raw,header=T,sep="\t",row.names=NULL)
gene_ids <- gex_table[,1]
cell_ids <- colnames(gex_table)[-1]
gex_table <- gex_table[,-1]

fetal_heart_test <- list(gex = list(NULL), pixel = list(NULL))

for (image_id in c(1:length(test_IDS))){
  
  coordinates <- gsub(meta_table[,1],pattern = "X",replacement = "")
  coordinate_file <- read.delim2(paste(main_path,"/images/",coordinate_list[image_id],sep=""),row.names=NULL)
  coordinate_file$coord <- paste(image_id,coordinate_file$x,coordinate_file$y,sep="x")
  
  coordinate_file <- coordinate_file[match(coordinates[grep(pattern = paste("^",image_id,"x",sep=""),coordinates)],coordinate_file$coord),]
  
  histo <- magick::image_read(paste(main_path,"/images/",hires_image[test_IDS[image_id]],sep=""))
  pixel_extracted <- extract_pixels(image = histo, coords = apply(coordinate_file[,c(5,6)],2,as.numeric))  
  fetal_heart_test$pixel[[image_id]] <- pixel_extracted
  fetal_heart_test$gex[[image_id]] <- t(gex_table[,paste("X",coordinates[grep(pattern = paste("^",image_id,"x",sep=""),coordinates)],sep="")])

}


save(fetal_heart_test, file = "./data/workflow/fetal_heart/fetal_heart_test.RData")

rm(fetal_heart_test)
gc()



fetal_heart_train <- list(gex = list(NULL), pixel = list(NULL))

for (image_id in c(1:length(train_IDS))){
  
  coordinates <- gsub(meta_table[,1],pattern = "X",replacement = "")
  coordinate_file <- read.delim2(paste(main_path,"/images/",coordinate_list[image_id],sep=""),row.names=NULL)
  coordinate_file$coord <- paste(image_id,coordinate_file$x,coordinate_file$y,sep="x")
  
  coordinate_file <- coordinate_file[match(coordinates[grep(pattern = paste("^",image_id,"x",sep=""),coordinates)],coordinate_file$coord),]
  
  histo <- magick::image_read(paste(main_path,"/images/",hires_image[train_IDS[image_id]],sep=""))
  pixel_extracted <- extract_pixels(image = histo, coords = apply(coordinate_file[,c(5,6)],2,as.numeric))  
  fetal_heart_train$pixel[[image_id]] <- pixel_extracted
  fetal_heart_train$gex[[image_id]] <- t(gex_table[,paste("X",coordinates[grep(pattern = paste("^",image_id,"x",sep=""),coordinates)],sep="")])
  
}


save(fetal_heart_train, file = "./data/workflow/fetal_heart/fetal_heart_train.RData")

gene_consensus <- gene_ids
save(gene_consensus,file = "./data/workflow/fetal_heart/gene_consensus.RData")

