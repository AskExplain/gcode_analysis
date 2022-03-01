load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/adipose_test.RData")

load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/gcode___adipose.5.all.models.Rdata")

library(lsa)

fa0 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$beta[[1]]
fb0 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$beta[[2]]
fi0 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$intercept[[2]]
fi1 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$intercept[[1]]
fc <- gcode.all.models$gcode.non_tumour[[1]]$main.code$code[[1]]

main_metrics <- c()
predicted_gex.alpha <- as.matrix(do.call('rbind',(ADIPOSE_test$pixel)) - c(fi0))%*%fb0%*%t(fc)%*%MASS::ginv(fc%*%t(fc))
predicted_gex <- predicted_gex.alpha%*%fc%*%t(fa0) +c(fi1)

observed_gex <- Matrix::t(do.call('cbind',(ADIPOSE_test$gex)))
row.names(predicted_gex) <- row.names(observed_gex)

gene_main_metrics <- rbind(main_metrics,do.call('rbind',pbmcapply::pbmclapply(c(1:dim(observed_gex)[2]),function(id){
  
  main_cor <- c(cosine(as.numeric(observed_gex[,id]),as.numeric(predicted_gex[,id])))
  main_mae <- c(mean(abs(as.numeric(observed_gex[,id])-as.numeric(predicted_gex[,id]))))
  
  list(main_cor = main_cor,main_mae = main_mae)
},mc.cores = 8)))

samples_main_metrics <- rbind(main_metrics,do.call('rbind',pbmcapply::pbmclapply(c(1:dim(observed_gex)[1]),function(id){
  
  main_cor <- c(cosine(as.numeric(observed_gex[id,]),as.numeric(predicted_gex[id,])))
  main_mae <- c(mean(abs(as.numeric(observed_gex[id,])-as.numeric(predicted_gex[id,]))))
  
  list(main_cor = main_cor,main_mae = main_mae)
},mc.cores = 8)))



if (plotting){

  X=1  
  par(mfcol=c(2,2))
  boxplot(unlist(samples_main_metrics[,1]),ylim=c(-1,1),main=median(unlist(samples_main_metrics[,1]),na.rm = T),xlab=X)
  boxplot(unlist(samples_main_metrics[,2]),main=median(unlist(samples_main_metrics[,2]),na.rm = T),xlab=X)
  
  boxplot(unlist(gene_main_metrics[,1]),ylim=c(-1,1),main=median(unlist(gene_main_metrics[,1]),na.rm = T),xlab=X)
  boxplot(unlist(gene_main_metrics[,2]),main=median(unlist(gene_main_metrics[,2]),na.rm = T),xlab=X)
  
  
  
  print(quantile(unlist(samples_main_metrics[,1]),c(0:20)/20))
  print(quantile(unlist(gene_main_metrics[,1]),c(0:20)/20,na.rm = T))
  
  
  
  
  for (i in order(unlist(gene_main_metrics[,1]),decreasing = T)[1:10]){
    plot(observed_gex[,i],predicted_gex[,i])
  }
  
  
  for (i in order(unlist(samples_main_metrics[,1]),decreasing = T)[1:10]){
    plot(observed_gex[i,],predicted_gex[i,])
  }
}