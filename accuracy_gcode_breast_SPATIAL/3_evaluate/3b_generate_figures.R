setwd("~/Documents/main_files/AskExplain/generative_encoder/")


plotting=FALSE
source("./github/accuracy_gcode_breast_SPATIAL/3_evaluate/3a_evaluate_gcode.R")
source("./github/accuracy_gcode_breast_SPATIAL/3_evaluate/plot_similarity_predict_vs_observe.R")

plot_similarity_predict_vs_observe(sample_similarity = samples_main_metrics, feature_similarity = gene_main_metrics, title_name = "", tissue_name = "breast")
