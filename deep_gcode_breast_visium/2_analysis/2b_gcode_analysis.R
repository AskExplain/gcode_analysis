###
# The R files for Generative Encoding can be found at
# https://drive.google.com/drive/folders/1a6txh_ZBcpPnLFj9mTknfwZ3A0dNTIlS
###

devtools::document("~/Documents/main_files/AskExplain/generative_encoder/package/generative_encoder_v2022.4/")
remove.packages("gcode")
devtools::install_local("~/Documents/main_files/AskExplain/generative_encoder/package/generative_encoder_v2022.4/")


load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/breast_visium_GE_data.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/BC_data_pixel.RData")



data_list <- list(ST_breast_gex_non.tumour = as.matrix(BC_data$non_tumour$gex),
                  ST_breast_pixel_non.tumour = as.matrix(BC_data_pixel)
)

rm(BC_data)
rm(BC_data_pixel)
gc()


for (main_dim in c(5,10,50,100,200,300,400,500)){
    
  config <- gcode::extract_config(F)
  config$init <- c("irlba","irlba")
  config$i_dim <- main_dim
  config$j_dim <- main_dim
  join <- list(alpha = c(1,1), beta = c(1,2), code = c(1,1) )
  config$bootstrap <- 2000 / main_dim
  config$tol <- 1
  config$regularise$a <- 0
  config$regularise$l <- 0.5
  
  gcode.non_tumour <- gcode::dgcode(data_list = data_list, config = config, join = join)
  
  gcode.all.models <- list(gcode.non_tumour = gcode.non_tumour)
  
  save(gcode.all.models,file=paste("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gcode___",main_dim,".all.models.Rdata",sep=""))
  
}
