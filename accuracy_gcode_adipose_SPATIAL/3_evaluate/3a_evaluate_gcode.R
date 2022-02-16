load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/gene_consensus.RData")
load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/adipose_test.RData")

load("~/Documents/main_files/AskExplain/generative_encoder/data/workflow/adipose/gcode___all.models.Rdata")

a0 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$beta[[1]]
b0 <- gcode.all.models$gcode.non_tumour[[1]]$main.parameters$beta[[2]]

main_metrics <- c()
predicted_gex <- as.matrix(do.call('rbind',ADIPOSE_test$pixel))%*%b0%*%t(sign(median(a0[a0!=0])*median(b0[b0!=0]))*a0)

observed_gex <- Matrix::t(do.call('cbind',(ADIPOSE_test$gex)))
row.names(predicted_gex) <- row.names(observed_gex)

gene_main_metrics <- rbind(main_metrics,do.call('rbind',pbmcapply::pbmclapply(c(1:dim(observed_gex)[2]),function(id){
  
  main_cor <- c(cor(as.numeric(observed_gex[,id]),as.numeric(predicted_gex[,id])))
  main_mae <- c(mean(abs(as.numeric(observed_gex[,id])-as.numeric(predicted_gex[,id]))))
  
  list(main_cor = main_cor,main_mae = main_mae)
},mc.cores = 8)))

samples_main_metrics <- rbind(main_metrics,do.call('rbind',pbmcapply::pbmclapply(c(1:dim(observed_gex)[1]),function(id){
  
  main_cor <- c(cor(as.numeric(observed_gex[id,]),as.numeric(predicted_gex[id,])))
  main_mae <- c(mean(abs(as.numeric(observed_gex[id,])-as.numeric(predicted_gex[id,]))))
  
  list(main_cor = main_cor,main_mae = main_mae)
},mc.cores = 8)))



if (plotting){
  
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