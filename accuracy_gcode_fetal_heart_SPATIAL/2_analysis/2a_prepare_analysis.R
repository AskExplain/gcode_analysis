setwd("~/Documents/main_files/AskExplain/generative_encoder/")
load("./data/workflow/fetal_heart/fetal_heart_train.RData")
load("./data/workflow/fetal_heart/fetal_heart_test.RData")

fetal_heart_data_pixel <- do.call('rbind',fetal_heart_train$pixel)

save(fetal_heart_data_pixel,file = "~/Documents/main_files/AskExplain/generative_encoder/data/workflow/fetal_heart/fetal_heart_data_pixel.RData")

rm(fetal_heart_data_pixel)
gc()

fetal_heart_data <- list(gex=NULL)
fetal_heart_data$gex <- do.call('rbind',fetal_heart_train$gex)

save(fetal_heart_data,file = "~/Documents/main_files/AskExplain/generative_encoder/data/workflow/fetal_heart/fetal_heart_visium_GE_data.RData")
