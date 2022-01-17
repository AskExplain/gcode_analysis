setwd("")
source("./mouse_hippocampus/1_prepare/load_hippocampus.R")

expression <- log(1+t(expression/(1e-3+colSums(expression)))) # log library normalisation
expression <- expression * 10e6 # rescaling 

library(gcode)

gcode.config <- gcode::extract_config(F)
gcode.config$method <- "gcode"
gcode.config$init <- c("irlba","irlba")
gcode.config$i_dim <- 100
gcode.config$j_dim <- 100
gcode.config$min_iter <- 3
gcode.config$max_iter <- 100
gcode.config$tol <- 1e-2

gcode.model <- gcode::gcode(
  data_list = list(as.matrix(expression)),
  config = gcode.config,
  join = list(alpha=1,beta=1,code=1)
)

v <- gcode.model$main.parameters$beta[[1]]
umap_ev <- umap::umap(as.matrix(expression)%*%(v))$layout
umap_e <- umap::umap(expression)$layout

pdf("./paper/method/2_adult_mouse_hippocampus_clusters.pdf",width = 15,height = 10)
par(mfrow=c(2,3))
plot(umap_ev,xlab="UMAP 1", ylab = "UMAP 2", main = "Adult Mouse Hippocampus via [gcode + UMAP]", col = as.factor(cell_clusters$group), pch = 19, cex = 0.8)
plot(umap_e,xlab="UMAP 1", ylab = "UMAP 2", main = "Adult Mouse Hippocampus via [UMAP]", col = as.factor(cell_clusters$group), pch = 19, cex = 0.8)
plot(0,0,col="white",xaxt="n",yaxt="n",ylab="",xlab="",type="n",bty="n")
legend("center", cex=2 ,legend = unique(cell_clusters$group), fill = as.factor(unique(cell_clusters$group)))

plot(umap_ev,xlab="UMAP 1", ylab = "UMAP 2", main = "Adult Mouse Hippocampus via [gcode + UMAP]", col = c(rev(colorRamps::primary.colors(length(unique(cell_clusters$group.1))))[as.integer(as.factor(cell_clusters$group.1))]), pch = 19, cex = 0.8)
plot(umap_e,xlab="UMAP 1", ylab = "UMAP 2", main = "Adult Mouse Hippocampus via [UMAP]", col = c(rev(colorRamps::primary.colors(length(unique(cell_clusters$group.1))))[as.integer(as.factor(cell_clusters$group.1))]), pch = 19, cex = 0.8)
plot(0,0,col="white",xaxt="n",yaxt="n",ylab="",xlab="",type="n",bty="n")
legend("center", cex=0.95 ,legend = unique(cell_clusters$group.1), fill = rev(colorRamps::primary.colors(length(unique(cell_clusters$group.1)))))
dev.off()

