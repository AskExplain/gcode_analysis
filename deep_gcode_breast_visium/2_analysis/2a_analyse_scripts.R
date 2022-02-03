
setwd("~/Documents/main_files/AskExplain/generative_encoder/")
load("./data/workflow/breast_visium.RData")

gene_consensus <- lapply(BC_all,function(X){
  colnames(X$gex$non_tumour)
})


gene_consensus <- Reduce("intersect",lapply(which(do.call('c',lapply(gene_consensus,function(X){
  length(X)!=0
}))),function(Y){
  gene_consensus[[Y]]
}))


BC_data_pixel <- do.call('rbind',lapply(BC_all,function(X){
  if (!is.null(X$pixel$non_tumour)){
    gcode::prepare_data((X$pixel$non_tumour),log = F,center = T,scale.norm = F)
  }
}))

save(BC_data_pixel,file = "~/Documents/main_files/AskExplain/generative_encoder/data/workflow/BC_data_pixel.RData")

rm(BC_data_pixel)
gc()


BC_data <- list(non_tumour=list(covariates=NULL,gex=NULL))

BC_data$non_tumour$gex <- do.call('rbind',lapply(c(BC_all),function(X){
  if (!is.null(X$gex$non_tumour)){
    gcode::prepare_data(X$gex$non_tumour[,gene_consensus],log = F,center = T,scale.norm = F)
  }
}))

BC_data$non_tumour$covariates <- do.call('rbind',lapply(c(BC_all),function(X){
  X$covariates$non_tumour
}))


rm(BC_all)
gc()


library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

IDs <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
             filters = "ensembl_gene_id", values = gene_consensus,
             mart = mart)

gene_consensus <- IDs

BC_data$non_tumour$gex <- BC_data$non_tumour$gex[,gene_consensus$ensembl_gene_id]

save(gene_consensus,file = "~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gene_consensus.RData")
save(BC_data,file = "~/Documents/main_files/AskExplain/generative_encoder/data/workflow/breast_visium_GE_data.RData")

