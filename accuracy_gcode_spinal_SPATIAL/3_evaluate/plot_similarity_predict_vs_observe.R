plot_similarity_predict_vs_observe <- function(sample_similarity,
                                               feature_similarity,
                                               title_name,
                                               tissue_name){
  
  library(ggplot2)
  
  lm <- rbind(c(1,2),
              c(1,2),
              c(1,2))
  
  g1 <- ggplot(data.frame(Measure=do.call('c',sample_similarity[,1]),Metric="Predicted vs Observed \n Gene-wise Pearson Correlation"), aes(x=Metric,y=Measure)) + 
    geom_violin() + ylim(-1,1) 
  
  g2 <- ggplot(data.frame(Measure=do.call('c',feature_similarity[,1]),Metric="Predicted vs Observed \n Sample-wise Pearson Correlation"), aes(x=Metric,y=Measure)) + 
    geom_violin() + ylim(-1,1) 
  
  gg_plots <- list(g1,g2)
  
  library(grid)
  library(gridExtra)
  final_plots <- arrangeGrob(
    grobs = gg_plots,
    layout_matrix = lm
  )
  
  
  ggsave(final_plots,filename = paste("./github/github_figures/jpeg_accuracy_gcode_",tissue_name,"_SPATIAL/",title_name,"accuracy_",tissue_name,"_spatial.png",sep=""),width = 7,height=11)
  
}