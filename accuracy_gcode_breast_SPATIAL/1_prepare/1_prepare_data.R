
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


setwd("~/Documents/main_files/AskExplain/generative_encoder/")
ids <- list.files("./data/external/breast/")
breast_cancer_ids <- substr(do.call('c',lapply(strsplit(ids[grep("Coords",ids)],split = "_"),function(X){paste(X[1],X[2],sep="_")})),3,10)

extract_pixels <- function(image,coords){
  
  coord_pixels <- do.call('rbind',pbmcapply::pbmclapply(c(1:dim(coords)[1]),function(i){
    X <- coords[i,1] - 50
    Y <- coords[i,2] - 50
    
    cropped_image <- magick::image_crop(image = image, geometry = paste("100x100+",X,"+",Y,sep=""))
    coord_pixels <- c(as.numeric(magick::image_data(cropped_image, 'rgb')))
    
  },mc.cores = 8))
  
  return(coord_pixels) 
}


BC_gex <- list(NULL)
BC_pixel <- list(NULL)

for (id in breast_cancer_ids){
  
  print(id)

  BC_covar <- list()
  
  all_files <- ids[grep(id,ids)]
  
  coords <- read.table(paste("./data/external/breast/",all_files[1],sep=""))
  gex <- read.table(paste("./data/external/breast/",all_files[2],sep=""))
  histo <- magick::image_read(paste("./data/external/breast/",all_files[3],sep=""))
  spots <- read.csv(paste("./data/external/breast/",all_files[4],sep=""))
  
  spots_internal <- cbind(spots[spots[,1]%in%row.names(gex),],id)
  BC_covar$non_tumour <- spots_internal[coords$tumor=="non",]
  
  if (length(BC_covar$non_tumour$X)>0){
    BC_gex[[id]] <- round(gex[row.names(gex)%in%BC_covar$non_tumour[,1],],3)
    BC_pixel[[id]] <- round(extract_pixels(image = histo, coords = BC_covar$non_tumour[,c(2,3)]),3)
  }
  
}

BC_all <- list(gex=BC_gex,pixel=BC_pixel)

save(BC_all, file = "./data/workflow/breast.RData")

