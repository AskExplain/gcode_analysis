setwd("~/Documents/main_files/AskExplain/generative_encoder/")

hires_image <- list.files("./data/external/spinal/images/HE/",recursive = T,pattern = "*.jpg")
coordinate_list <- list.files("./data/external/spinal/annotations/",recursive = T,pattern = "*.tsv")
gex_raw <- list.files("./data/external/spinal/count_matrices/",recursive = T,pattern = "*.txt")

image_IDS <- do.call('c',lapply(strsplit(coordinate_list,split = "\\."),function(X){X[1]}))

main_path <- c("~/Documents/main_files/AskExplain/generative_encoder/data/external/spinal/")

extract_pixels <- function(image,coords){
  
  coord_pixels <- do.call('rbind',pbmcapply::pbmclapply(c(1:dim(coords)[1]),function(i){
    X <- coords[i,1] - 32
    Y <- coords[i,2] - 32
    
    cropped_image <- magick::image_crop(image = image, geometry = paste("64x64+",X,"+",Y,sep=""))
    coord_pixels <- c(as.numeric(magick::image_data(cropped_image, 'rgb')))
    
  },mc.cores = 8))
  
  return(coord_pixels) 
}

set.seed(1)
train_IDS <- sample(c(1:length(image_IDS)),length(image_IDS)*0.6)
test_IDS <- c(1:length(image_IDS))[-train_IDS]


spinal_test <- list(gex = list(NULL), pixel = list(NULL))


for (image_id in c(1:length(test_IDS))){
  
  spinal_test$gex[[image_id]] <- t(read.delim2(paste(main_path,"/count_matrices/",gex_raw[test_IDS[image_id]],sep=""),header=T,sep="\t",row.names=NULL))
  colnames(spinal_test$gex[[image_id]]) <- spinal_test$gex[[image_id]][1,]
  row_ids <- row.names(spinal_test$gex[[image_id]])[-1]
  spinal_test$gex[[image_id]] <- apply(spinal_test$gex[[image_id]][-1,],2,as.numeric)
  gene_ids <- colnames(spinal_test$gex[[image_id]])
  
  coordinate_file <- apply(do.call('rbind',strsplit(gsub(row_ids,pattern = "X",replacement = ""),split = "_")),2,as.numeric)
  
  histo <- magick::image_read(paste(main_path,"/images/HE/",hires_image[test_IDS[image_id]],sep=""))
  main_coords <- apply(as.matrix(apply(coordinate_file,2,as.numeric))*194/6400,2,function(X){pmin(X,1)})*as.numeric(magick::image_info(histo)[2:3])
  pixel_extracted <- extract_pixels(image = histo, coords = main_coords)  
  
  if (!is.character(pixel_extracted)){
    spinal_test$pixel[[image_id]] <- pixel_extracted
  } else {
    spinal_test$pixel[[image_id]] <- NULL
    spinal_test$gex[[image_id]] <- NULL
  }
  
}


save(spinal_test, file = "./data/workflow/spinal/spinal_test.RData")

rm(spinal_test)
gc()



spinal_train <- list(gex = list(NULL), pixel = list(NULL))

for (image_id in c(1:length(train_IDS))){
  
  spinal_train$gex[[image_id]] <- t(read.delim2(paste(main_path,"/count_matrices/",gex_raw[train_IDS[image_id]],sep=""),header=T,sep="\t",row.names=NULL))
  colnames(spinal_train$gex[[image_id]]) <- spinal_train$gex[[image_id]][1,]
  row_ids <- row.names(spinal_train$gex[[image_id]])[-1]
  
  spinal_train$gex[[image_id]] <- apply(spinal_train$gex[[image_id]][-1,],2,as.numeric)
  gene_ids <- colnames(spinal_train$gex[[image_id]])
  
  coordinate_file <- apply(do.call('rbind',strsplit(gsub(row_ids,pattern = "X",replacement = ""),split = "_")),2,as.numeric)
  
  histo <- magick::image_read(paste(main_path,"/images/HE/",hires_image[train_IDS[image_id]],sep=""))
  main_coords <- apply(as.matrix(apply(coordinate_file,2,as.numeric))*194/6400,2,function(X){pmin(X,1)})*as.numeric(magick::image_info(histo)[2:3])
  pixel_extracted <- extract_pixels(image = histo, coords = main_coords)  
  
  if (!is.character(pixel_extracted)){
    spinal_train$pixel[[image_id]] <- pixel_extracted
  } else {
    spinal_train$pixel[[image_id]] <- NULL
    spinal_train$gex[[image_id]] <- NULL
  }
  
}


save(spinal_train, file = "./data/workflow/spinal/spinal_train.RData")
