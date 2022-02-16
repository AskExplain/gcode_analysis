###
# The R files for Generative Encoding can be found at
# https://drive.google.com/drive/folders/1a6txh_ZBcpPnLFj9mTknfwZ3A0dNTIlS
###

library(gcode)

setwd("~/Documents/main_files/AskExplain/generative_encoder")

load("./data/workflow/gene_consensus.RData")
load("./data/workflow/breast_visium_GE_data.RData")
load("./data/workflow/BC_data_pixel.RData")

row.names(BC_data_pixel) <- row.names(BC_data$non_tumour$gex)

set.seed(1)
train_IDS <- sample(c(1:dim(BC_data_pixel)[1]),dim(BC_data_pixel)[1]*0.70)

train_data_list <- list(GEX=as.matrix(BC_data$non_tumour$gex)[train_IDS,],PIXEL_TRAIN=BC_data_pixel[train_IDS,])
test_data_list <- list(GEX=as.matrix(BC_data$non_tumour$gex)[-train_IDS,],PIXEL_TRAIN=BC_data_pixel[-train_IDS,])

save(train_data_list,file = "./data/workflow/BREAST_train.RData")
save(test_data_list,file = "./data/workflow/BREAST_test.RData")

load("./data/workflow/BREAST_train.RData")

config <- gcode::extract_config(F)
config$init <- c("irlba","irlba")
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

gcode.non_tumour <- gcode::dgcode(data_list = train_data_list, config = config, join = join)

gcode.all.models <- list(gcode.non_tumour=gcode.non_tumour)

save(gcode.all.models,file=paste("./data/workflow/gcode___all.models.Rdata",sep=""))

