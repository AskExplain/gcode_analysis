setwd("~/Documents/main_files/AskExplain/generative_encoder/")

load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/spinal/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/spinal/spinal_visium_GE_data.RData")


data_list <- spinal_data


for (dim_all in c(5)){
  config <- gcode::extract_config(F)
  config$init <- c("irlba","irlba")
  config$transform$log <- F
  config$transform$center <- F
  config$transform$norm <- F
  config$i_dim <- dim_all
  config$j_dim <- dim_all
  config$tol <- 1
  config$regularise$a <- 0
  config$regularise$l <- 0
  config$dimension_reduction <- F
  
  
  #############
  # IMPORTANT #
  #############
  
  # This join$complete list can get complicated
  # The reason for this is to model different ASPECTS of the data (technical platform biases, samples that share the same underlying biological functional mechanism etc.)
  # join$complete$data_list - two datasets in the main data list, each identifier here represents a single dataset
  # join$complete$alpha - represents both technical biases and common signals between same samples of different modalities (
  #     alpha = 1 for platform specific biases in gene expression measurement, and,
  #     alpha = 2 for observed sample space of gene expression and histology, and, 
  #     alpha = 3 for platform bias correction on observed sample space of gene expression and histology
  #     alpha = 4 for platform specific biases in measuring spatial histology
  # join$complete$beta - represents the true gene expression and histology signals (
  #     beta = 1 for true gene expression signal, and,
  #     beta = 2 for true spatial histology signal
  # join$complete$code - represents the commonality between datasets (
  #     code = 1 assigns commonality to both platform specific biases of measurements for gene expression and observations (for correction in code = 2)
  #     code = 2 corrects observations for confounding biases
  #     code = 3 assigns commonality to both platform specific biases of measurements for spatial histology and observations (for correction in code = 2)
  join <- list(complete = list(data_list = c(1,2),
                               alpha = c(1,1),
                               beta = c(1,2),
                               code = c(1,1)
  ),
  labels = list(alpha=NULL,
                beta=NULL))
  
  references <- gcode::extract_references_framework(F)
  references$data_list <- c(1,0)
  
  gcode.non_tumour <- gcode::gcode(data_list = data_list, config = config, join = join, references = references)
  
  gcode.all.models <- list(gcode.non_tumour=list(gcode.non_tumour))
  
  save(gcode.all.models,file=paste(paste("./data/workflow/spinal/gcode___spinal.",dim_all,".all.models.Rdata",sep="")))
}
