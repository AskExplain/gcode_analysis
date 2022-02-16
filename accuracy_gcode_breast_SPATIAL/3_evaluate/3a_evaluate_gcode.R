setwd("~/Documents/main_files/AskExplain/generative_encoder/")

load("./data/workflow/gene_consensus.RData")
load("./data/workflow/BREAST_test.RData")

load("./data/workflow/gcode___all.models.Rdata")

a0 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$beta[[1]]
b0 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$beta[[2]]

main_metrics <- c()
predicted_gex <- test_data_list$PIXEL_TRAIN%*%b0%*%t(sign(median(a0[a0!=0])*median(b0[b0!=0]))*a0)

observed_gex <- (test_data_list$GEX)
row.names(predicted_gex) <- row.names(observed_gex)

gene_main_metrics <- rbind(main_metrics,do.call('rbind',pbmcapply::pbmclapply(c(1:dim(observed_gex)[2]),function(id){
  
  ids_keep <- which(observed_gex[,id]!=0)
  main_cor <- c(cor(as.numeric(observed_gex[ids_keep,id]),as.numeric(predicted_gex[ids_keep,id])))
  main_mae <- c(mean(abs(as.numeric(observed_gex[ids_keep,id])-as.numeric(predicted_gex[ids_keep,id]))))
  
  list(main_cor = main_cor,main_mae = main_mae)
},mc.cores = 8)))

samples_main_metrics <- rbind(main_metrics,do.call('rbind',pbmcapply::pbmclapply(c(1:dim(observed_gex)[1]),function(id){
  
  ids_keep <- which(observed_gex[id,]!=0)
  main_cor <- c(cor(as.numeric(observed_gex[id,ids_keep]),as.numeric(predicted_gex[id,ids_keep])))
  main_mae <- c(mean(abs(as.numeric(observed_gex[id,ids_keep])-as.numeric(predicted_gex[id,ids_keep]))))
  
  list(main_cor = main_cor,main_mae = main_mae)
},mc.cores = 8)))



if (plotting){
  par(mfcol=c(2,2))
  boxplot(unlist(samples_main_metrics[,1]),ylim=c(-1,1),main=median(unlist(samples_main_metrics[,1]),na.rm = T),xlab=X)
  boxplot(unlist(samples_main_metrics[,2]),main=median(unlist(samples_main_metrics[,2]),na.rm = T),xlab=X)
  
  boxplot(unlist(gene_main_metrics[,1]),ylim=c(-1,1),main=median(unlist(gene_main_metrics[,1]),na.rm = T),xlab=X)
  boxplot(unlist(gene_main_metrics[,2]),main=median(unlist(gene_main_metrics[,2]),na.rm = T),xlab=X)
  
  
  print(quantile(unlist(samples_main_metrics[,1]),c(0:20)/20,na.rm = T))
  print(quantile(unlist(gene_main_metrics[,1]),c(0:20)/20,na.rm = T))
  
  
  for (i in order(unlist(gene_main_metrics[,1]),decreasing = T)[1:10]){
    plot(observed_gex[,i],predicted_gex[,i])
  }
  
  
  for (i in order(unlist(samples_main_metrics[,1]),decreasing = T)[1:10]){
    plot(observed_gex[i,],predicted_gex[i,])
  }
  
}

