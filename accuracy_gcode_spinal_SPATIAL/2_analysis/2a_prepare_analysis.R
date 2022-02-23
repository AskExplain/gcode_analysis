setwd("~/Documents/main_files/AskExplain/generative_encoder/")
load("./data/workflow/spinal/spinal_train.RData")
load("./data/workflow/spinal/spinal_test.RData")

gene_ids <- lapply(c(spinal_train$gex,spinal_test$gex),colnames)
gene_ids <- lapply(which(!do.call('c',lapply(gene_ids,is.null))),function(X){gene_ids[[X]]})
gene_ids <- lapply(which(do.call('c',lapply(gene_ids,function(X){X[1]=="Gnai3"}))),function(X){gene_ids[[X]]})
gene_consensus <- Reduce("intersect",gene_ids)
save(gene_consensus,file = "./data/workflow/spinal/gene_consensus.RData")




spinal_data <- list(gex=NULL,pixel=NULL)
spinal_data$gex <- do.call('rbind',lapply(c(1:length(spinal_train$gex)),function(X){
  print(X)
  if (!is.null(spinal_train$gex[[X]])){
    if (colnames(spinal_train$gex[[X]])[1]=="Gnai3"){
    spinal_train$gex[[X]][,gene_consensus]
    }
  }
}))

spinal_data$pixel <- do.call('rbind',lapply(c(1:length(spinal_train$gex)),function(X){
  print(X)
  if (!is.null(spinal_train$gex[[X]])){
    if (colnames(spinal_train$gex[[X]])[1]=="Gnai3"){
      spinal_train$pixel[[X]]
    }
  }
}))

save(spinal_data,file = "~/Documents/main_files/AskExplain/generative_encoder/data/workflow/spinal/spinal_visium_GE_data.RData")

rm(spinal_data,spinal_train);gc()



spinal_data_test <- list(gex=NULL,pixel=NULL)
spinal_data_test$gex <- do.call('rbind',lapply(c(1:length(spinal_test$gex)),function(X){
  print(X)
  if (!is.null(spinal_test$gex[[X]])){
    if (colnames(spinal_test$gex[[X]])[1]=="Gnai3"){
      spinal_test$gex[[X]][,gene_consensus]
    }
  }
}))

spinal_data_test$pixel <- do.call('rbind',lapply(c(1:length(spinal_test$gex)),function(X){
  print(X)
  if (!is.null(spinal_test$gex[[X]])){
    if (colnames(spinal_test$gex[[X]])[1]=="Gnai3"){
      spinal_test$pixel[[X]]
    }
  }
}))

save(spinal_data_test,file = "~/Documents/main_files/AskExplain/generative_encoder/data/workflow/spinal/spinal_visium_GE_data_test.Rdata")
