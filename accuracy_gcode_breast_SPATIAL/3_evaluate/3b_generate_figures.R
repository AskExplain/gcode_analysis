setwd("~/Documents/main_files/AskExplain/generative_encoder/")

plotting=FALSE
source("./github/accuracy_gcode_breast_SPATIAL/3_evaluate/3a_evaluate_gcode.R")


library(ggplot2)

lm <- rbind(c(1,2,2,2),
            c(1,2,2,2),
            c(1,2,2,2))


gene_consensus$labels <- gene_consensus$hgnc_symbol
gene_consensus$labels[-order(predicted_gex[order(do.call('c',samples_main_metrics[,1]),decreasing = T)[1],],decreasing=T)[1:30]] <- c("")


g1 <- ggplot(data.frame(Measure=do.call('c',samples_main_metrics[,1]),Metric="Gene-wise Pearson Correlation"), aes(x=Metric,y=Measure)) + 
  geom_violin() + ylim(-1,1)

g2 <- ggplot(data.frame(Observed=observed_gex[order(do.call('c',samples_main_metrics[,1]),decreasing = T)[1],],
                        Predicted=predicted_gex[order(do.call('c',samples_main_metrics[,1]),decreasing = T)[1],])) + 
               aes(x=Observed, y=Predicted, label = gene_consensus$labels) + geom_point() + ggrepel::geom_label_repel()

gg_plots <- list(g1,g2)

library(grid)
library(gridExtra)
final_plots <- arrangeGrob(
  grobs = gg_plots,
  layout_matrix = lm
)


ggsave(final_plots,filename = "./github/github_figures/accuracy_breast_spatial.png",width = 15,height=10)

       