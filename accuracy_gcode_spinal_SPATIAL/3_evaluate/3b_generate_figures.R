setwd("~/Documents/main_files/AskExplain/generative_encoder/")

plotting=FALSE
source("./github/accuracy_gcode_spinal_SPATIAL/3_evaluate/3a_evaluate_gcode.R")

selected_subsample <- sample(c(1:dim(spinal_data_test$gex)[1]),size = 300)
cor_spinal_gene_wise <- cor(t(as.matrix(spinal_data_test$gex)[selected_subsample,]),t(as.matrix(spinal_data_test$gex)[-selected_subsample,][1:300,]))

selected_subsample <- sample(c(1:dim(spinal_data_test$gex)[2]),size = 300)
cor_spinal_sample_wise <- cor((as.matrix(spinal_data_test$gex)[,selected_subsample]),(as.matrix(spinal_data_test$gex)[,-selected_subsample][,1:300]))



library(ggplot2)

lm <- rbind(c(1,2,5,5,5,3,4),
            c(1,2,5,5,5,3,4),
            c(1,2,5,5,5,3,4))

gene_consensus <- data.frame(labels=colnames(predicted_gex))
gene_consensus$labels[-order(predicted_gex[order(do.call('c',samples_main_metrics[,1]),decreasing = T)[1],],decreasing=T)[1:100]] <- c("")


g1 <- ggplot(data.frame(Measure=do.call('c',samples_main_metrics[,1]),Metric="Predicted vs Observed \n Gene-wise Pearson Correlation"), aes(x=Metric,y=Measure)) + 
  geom_violin() + ylim(-1,1) + 
  ggtitle(label = paste("P-value: ",t.test(y=c(cor_spinal_gene_wise),x=do.call('c',samples_main_metrics[,1]))$p.value,sep=""))

g2 <- ggplot(data.frame(Measure=c(cor_spinal_gene_wise),Metric="Observed (random set 1) vs Observed (random set 2) \n Gene-wise Pearson Correlation"), aes(x=Metric,y=Measure)) + 
  geom_violin() + ylim(-1,1) +
  ggtitle(label = paste("Test-statistic: ",t.test(y=c(cor_spinal_gene_wise),x=do.call('c',samples_main_metrics[,1]))$statistic,sep=""))


g3 <- ggplot(data.frame(Measure=do.call('c',gene_main_metrics[,1]),Metric="Predicted vs Observed \n Sample-wise Pearson Correlation"), aes(x=Metric,y=Measure)) + 
  geom_violin() + ylim(-1,1) + 
  ggtitle(label = paste("P-value: ",t.test(y=c(cor_spinal_sample_wise),x=do.call('c',gene_main_metrics[,1]))$p.value,sep=""))


g4 <- ggplot(data.frame(Measure=c(cor_spinal_sample_wise),Metric="Observed (random set 1) vs Observed (random set 2) \n Sample-wise Pearson Correlation"), aes(x=Metric,y=Measure)) + 
  geom_violin() + ylim(-1,1) + 
  ggtitle(label = paste("Test-statistic: ",t.test(y=c(cor_spinal_sample_wise),x=do.call('c',gene_main_metrics[,1]))$statistic,sep=""))


g5 <- ggplot(data.frame(Observed=observed_gex[order(do.call('c',samples_main_metrics[,1]),decreasing = T)[1],],
                        Predicted=predicted_gex[order(do.call('c',samples_main_metrics[,1]),decreasing = T)[1],])) + 
  aes(x=Observed, y=Predicted, label = gene_consensus$labels) + geom_point() + ggrepel::geom_label_repel()


gg_plots <- list(g1,g2,g3,g4,g5)

library(grid)
library(gridExtra)
final_plots <- arrangeGrob(
  grobs = gg_plots,
  layout_matrix = lm
)


ggsave(final_plots,filename = "./github/github_figures/jpeg_accuracy_gcode_spinal_SPATIAL/accuracy_spinal_spatial.png",width = 24,height=10)

