
######## 	 All files can be found below:

### Spatial Transcriptomics dataset
# Data downloaded at 
# https://data.mendeley.com/datasets/29ntw7sh4r/5 and put in a file at ./data/external/breast_visium/
###

###
# The R files for Generative Encoding can be found at
# https://drive.google.com/drive/folders/1a6txh_ZBcpPnLFj9mTknfwZ3A0dNTIlS
###

#########


setwd("~/Documents/main_files/AskExplain/generative_encoder/")
ids <- list.files("./data/external/breast_visium/")
breast_cancer_ids <- substr(do.call('c',lapply(strsplit(ids[grep("Coords",ids)],split = "_"),function(X){paste(X[1],X[2],sep="_")})),3,10)


set.seed(100)

ids_train___ <- sample(c(1:length(breast_cancer_ids)),ceiling(0.9*length(breast_cancer_ids)))
train___ <- breast_cancer_ids[ids_train___]
test___ <- breast_cancer_ids[-ids_train___]

extract_pixels <- function(image,coords){
  
  coord_pixels <- c()
  
  for (i in 1:dim(coords)[1]){
    
    X <- coords[i,1] - 32
    Y <- coords[i,2] - 32
    
    cropped_image <- magick::image_crop(image = image, geometry = paste("64x64+",X,"+",Y,sep=""))
    coord_pixels <- rbind(coord_pixels,c(as.numeric(magick::image_data(cropped_image, 'rgb'))))
  }
 
  return(coord_pixels) 
}
  


BC_all <- list()

for (id in train___){
  
  print(id)

  BC_gex <- list()
  BC_pixel <- list()
  BC_covar <- list()
  
  all_files <- ids[grep(id,ids)]
  
  coords <- read.table(paste("./data/external/breast_visium/",all_files[1],sep=""))
  gex <- read.table(paste("./data/external/breast_visium/",all_files[2],sep=""))
  histo <- magick::image_read(paste("./data/external/breast_visium/",all_files[3],sep=""))
  spots <- read.csv(paste("./data/external/breast_visium/",all_files[4],sep=""))
  
  spots_internal <- cbind(spots[spots[,1]%in%row.names(gex),],id)
  BC_covar$tumour <- spots_internal[coords$tumor=="tumor",]
  BC_covar$border <- if (dim(spots_internal[coords$tumor=="non",])[1]>0){
    do.call('rbind',lapply(c(1:dim(spots_internal[coords$tumor=="non",])[1]),function(X){
    if (any(apply(c(spots_internal[coords$tumor=="non",][X,][,c(2,3)]) - BC_covar$tumour[,c(2,3)],1,function(Y){sqrt(sum((Y)^2))})<350)){
      return(spots_internal[coords$tumor=="non",][X,])
    } else {
      NULL
    }
    }))
    } else {
    NULL
    }
  BC_covar$non_tumour <- if (dim(spots_internal[coords$tumor=="non",])[1]>0){
    do.call('rbind',lapply(c(1:dim(spots_internal[coords$tumor=="non",])[1]),function(X){
    if (!any(apply(c(spots_internal[coords$tumor=="non",][X,][,c(2,3)]) - BC_covar$tumour[,c(2,3)],1,function(Y){sqrt(sum((Y)^2))})<350)){
      return(spots_internal[coords$tumor=="non",][X,])
    } else {
      NULL
    }
  }))
  } else {
    NULL
  }
  
  # if (length(BC_covar$tumour)>0){
  #   BC_gex$tumour <- gex[row.names(gex)%in%BC_covar$tumour[,1],]
  #   BC_pixel$tumour <- extract_pixels(image = histo, coords = BC_covar$tumour[,c(2,3)])
  # }
  # if (length(BC_covar$border)>0){
  #   BC_gex$border <- gex[row.names(gex)%in%BC_covar$border[,1],]
  #   BC_pixel$border <- extract_pixels(image = histo, coords = BC_covar$border[,c(2,3)])
  # }
  if (length(BC_covar$non_tumour)>0){
    BC_gex$non_tumour <- gex[row.names(gex)%in%BC_covar$non_tumour[,1],]
    BC_pixel$non_tumour <- extract_pixels(image = histo, coords = BC_covar$non_tumour[,c(2,3)])
  }
  
  BC_all[[id]] <- list(covariates = list(non_tumour = BC_covar$non_tumour), gex = list(non_tumour = BC_gex$non_tumour), pixel = list(non_tumour = BC_pixel$non_tumour))
  
}



save(BC_all, file = "./data/workflow/breast_visium.RData")
