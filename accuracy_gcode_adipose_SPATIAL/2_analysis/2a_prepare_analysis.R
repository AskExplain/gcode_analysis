setwd("~/Documents/main_files/AskExplain/generative_encoder/")
load("./data/workflow/adipose/adipose_train.RData")
load("./data/workflow/adipose/adipose_test.RData")

ADIPOSE_data_pixel <- do.call('rbind',ADIPOSE_train$pixel)

save(ADIPOSE_data_pixel,file = "~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/ADIPOSE_data_pixel.RData")

rm(ADIPOSE_data_pixel)
gc()

ADIPOSE_data <- list(gex=NULL)
ADIPOSE_data$gex <- t(do.call('cbind',ADIPOSE_train$gex))
ADIPOSE_data$gex <- ADIPOSE_data$gex / rowSums(ADIPOSE_data$gex)

save(ADIPOSE_data,file = "~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/adipose_visium_GE_data.RData")
