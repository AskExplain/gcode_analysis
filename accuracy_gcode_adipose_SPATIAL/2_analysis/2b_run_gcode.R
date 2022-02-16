###
# The R files for Generative Encoding can be found at
# https://drive.google.com/drive/folders/1a6txh_ZADIPOSEpPnLFj9mTknfwZ3A0dNTIlS
###

devtools::document("~/Documents/main_files/AskExplain/generative_encoder/package/generative_encoder_v2022.4/")
remove.packages("gcode")
devtools::install_local("~/Documents/main_files/AskExplain/generative_encoder/package/generative_encoder_v2022.4/")


load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/adipose_visium_GE_data.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/ADIPOSE_data_pixel.RData")

row.names(ADIPOSE_data_pixel) <- row.names(ADIPOSE_data$gex)

data_list <- list(GEX=as.matrix(ADIPOSE_data$gex),PIXEL_TRAIN=ADIPOSE_data_pixel)

rm(ADIPOSE_data)
rm(ADIPOSE_data_pixel)
gc()


config <- gcode::extract_config(F)
config$init <- list("irlba","irlba")
config$i_dim <- c(300)
config$j_dim <- 1000
config$regularise$a <- 0
config$regularise$l <- 0
config$dimension_reduction <- F

join <- list(complete = list(alpha = c(1,1),
                             beta = c(1,2),
                             code = c(1,1)
),
labels = list(alpha=NULL,
              beta=NULL))

gcode.non_tumour <- gcode::dgcode(data_list = data_list, config = config, join = join)

gcode.all.models <- list(gcode.non_tumour=gcode.non_tumour)

save(gcode.all.models,file=paste("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/gcode___all.models.Rdata",sep=""))


