setwd("~/Documents/main_files/AskExplain/generative_encoder/")

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
  
  coord_pixels <- do.call('rbind',pbmcapply::pbmclapply(c(1:dim(coords)[1]),function(i){
    X <- coords[i,1] - 50
    Y <- coords[i,2] - 50
    
    cropped_image <- magick::image_crop(image = image, geometry = paste("100x100+",X,"+",Y,sep=""))
    coord_pixels <- c(as.numeric(magick::image_data(cropped_image, 'rgb')))
    
  },mc.cores = 3))
  
  return(coord_pixels) 
}

set.seed(1)
train_IDS <- sample(c(1:length(image_IDS)),length(image_IDS)*0.6)
test_IDS <- c(1:length(image_IDS))[-train_IDS]


ADIPOSE_test <- list(gex = list(NULL), pixel = list(NULL))

for (image_id in c(1:length(test_IDS))){
  
  coordinate_file <- read.csv(paste(main_path,coordinate_list[test_IDS[image_id]],sep=""),header=F)
  scaling_file <- read.delim(paste(main_path,scaling_json[test_IDS[image_id]],sep=""),header=F)
  hires_scaling_factor <- as.numeric(gsub(" tissue_hires_scalef: ","",strsplit(scaling_file[1,1],split = ",")[[1]][2]))
  coordinate_file[,c(5,6)] <- coordinate_file[,c(5,6)]*hires_scaling_factor
  
  ADIPOSE_test$gex[[image_id]] <- (Matrix::readMM(paste(main_path,gex_raw[test_IDS[image_id]],sep="")))
  gene_ids <- read.table(paste(main_path,gex_features[test_IDS[image_id]],sep=""))
  row.names(ADIPOSE_test$gex[[image_id]]) <- gene_ids$V2
  
  IDS_remove <- which(Matrix::colSums(ADIPOSE_test$gex[[image_id]])<2000)
  
  ADIPOSE_test$gex[[image_id]] <- ADIPOSE_test$gex[[image_id]][,-IDS_remove]
  histo <- magick::image_read(paste(main_path,hires_image[test_IDS[image_id]],sep=""))
  ADIPOSE_test$pixel[[image_id]] <- extract_pixels(image = histo, coords = coordinate_file[-IDS_remove,c(5,6)])
  
}

save(ADIPOSE_test, file = "../workflow/adipose_test.RData")

rm(ADIPOSE_test)
gc()



ADIPOSE_train <- list(gex = list(NULL), pixel = list(NULL))

for (image_id in c(1:length(train_IDS))){
  
  coordinate_file <- read.csv(paste(main_path,coordinate_list[train_IDS[image_id]],sep=""),header=F)
  scaling_file <- read.delim(paste(main_path,scaling_json[train_IDS[image_id]],sep=""),header=F)
  hires_scaling_factor <- as.numeric(gsub(" tissue_hires_scalef: ","",strsplit(scaling_file[1,1],split = ",")[[1]][2]))
  coordinate_file[,c(5,6)] <- coordinate_file[,c(5,6)]*hires_scaling_factor
  
  ADIPOSE_train$gex[[image_id]] <- as.matrix(Matrix::readMM(paste(main_path,gex_raw[train_IDS[image_id]],sep="")))
  gene_ids <- read.table(paste(main_path,gex_features[train_IDS[image_id]],sep=""))
  row.names(ADIPOSE_train$gex[[image_id]]) <- gene_ids$V2
  
  IDS_remove <- which(Matrix::colSums(ADIPOSE_train$gex[[image_id]])<2000)
  
  ADIPOSE_train$gex[[image_id]] <- ADIPOSE_train$gex[[image_id]][,-IDS_remove]
  histo <- magick::image_read(paste(main_path,hires_image[train_IDS[image_id]],sep=""))
  ADIPOSE_train$pixel[[image_id]] <- extract_pixels(image = histo, coords = coordinate_file[-IDS_remove,c(5,6)])
  
}


save(ADIPOSE_train, file = "../workflow/adipose_train.RData")


gene_consensus <- gene_ids

save(gene_consensus,file = "../workflow/adipose/gene_consensus.RData")

