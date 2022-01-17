###
# The R files for Generative Encoding can be found at
# https://drive.google.com/drive/folders/1a6txh_ZBcpPnLFj9mTknfwZ3A0dNTIlS
###

devtools::document("~/Documents/main_files/AskExplain/generative_encoder/package/generative_encoder_v2022.1//")
remove.packages("gcode")
devtools::install_local("~/Documents/main_files/AskExplain/generative_encoder/package/generative_encoder_v2022.1/")


load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/breast_visium_GE_data.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/BC_data_pixel.RData")



data_list <- list(ST_breast_gex_non.tumour = as.matrix(BC_data$non_tumour$gex),
                  ST_breast_pixel_non.tumour = as.matrix(BC_data_pixel)
)

rm(BC_data)
rm(BC_data_pixel)
gc()



config <- gcode::extract_config(F)
config$init <- c("irlba","irlba")
config$i_dim <- 500
config$j_dim <- 500

join <- list(alpha = c(1,2), beta = c(1,2), code = c(1,1) )

gcode.non_tumour <- gcode::gcode(data_list = data_list, config = config, join = join)

gcode.all.models <- list(gcode.non_tumour = gcode.non_tumour)

save(gcode.all.models,file="~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gcode.all.models.Rdata")


