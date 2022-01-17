# Simulation sourced from
# https://rstudio-pubs-static.s3.amazonaws.com/222019_a2c617e90530422285420091c090d415.html


devtools::install_github("AskExplain/gcode")



# number of data points
k = 1000

# xf
set.seed(7)
x <- rnorm(n = k, mean = 0, sd = 1)

# random noise
set.seed(9)
e <- rnorm(n = k, mean = 0, sd = 1)

# y
y <- abs(x) + e

# Matrix A
# Store x and y into a 2 x k matrix
A <- matrix(c(x, y), byrow = TRUE, nrow = 2)

# matrix B
# Arbitrarily rotate, scale, and translate A 
theta <- pi / 4
rot.mtx <- matrix(c(sin(theta), -cos(theta), cos(theta), sin(theta)), ncol=2)
B <- 0.5*rot.mtx %*% A - 6


config <- gcode::extract_config(F)
config$init <- c("rsvd","rsvd")
config$i_dim <- 2
config$j_dim <- 2
A.B_cv.gcode <- gcode::gcode(data_list = list(as.matrix(t(A)),as.matrix(t(B))), join = list(alpha=c(1,2),beta = c(1,2), code=c(1,1)), config = config)



procrustes <- function(A, B){
  # center and normalize A 
  A.centered <- t(scale(t(A), center = TRUE, scale = FALSE))
  A.size <- norm(A.centered, type = "F") / (ncol(A) * nrow(A))
  A.normalized <- A.centered / A.size
  
  # center and normalize B
  B.centered <- t(scale(t(B), center = TRUE, scale = FALSE))
  B.size <- norm(B.centered, type = "F") / (ncol(B) * nrow(B))
  B.normalized <- B.centered / B.size
  
  # Rotation matrix T 
  svd.results <- svd(B.normalized %*% t(A.normalized))
  U <- svd.results$u
  V <- svd.results$v
  T <- V %*% t(U)
  
  # B transformed
  B.transformed <- T %*% B.normalized
  
  # Error after superimposition
  RSS <- norm(A.normalized - B.transformed,  type = "F")
  
  # Return
  return(list(A.normalized = A.normalized, B.normalized = B.normalized, rotation.mtx = T, B.transformed = B.transformed, RSS = RSS))
}

procrustes.results <- procrustes(A, B)

A.gcode <- as.data.frame(t(A))
B.gcode <- as.data.frame(t(A.B_cv.gcode$main.parameters$alpha[[1]])%*%(A.B_cv.gcode$main.parameters$alpha[[2]])%*%(t(B))%*%(A.B_cv.gcode$main.parameters$beta[[2]])%*%t(A.B_cv.gcode$main.parameters$beta[[1]]))

A.gcode_samples <- as.data.frame(t(A))
B.gcode_samples <- as.data.frame(t(A.B_cv.gcode$main.parameters$alpha[[1]])%*%(A.B_cv.gcode$main.parameters$alpha[[2]])%*%(t(B)))

A.gcode_features <- as.data.frame(t(A))
B.gcode_features <- as.data.frame((t(B))%*%(A.B_cv.gcode$main.parameters$beta[[2]])%*%t(A.B_cv.gcode$main.parameters$beta[[1]]))

A.pc.gcode <- (t(procrustes.results$A.normalized))
B.pc.gcode <- (t(procrustes.results$B.transformed))


data.gcode <- rbind(t(A), t(B), A.gcode_samples, B.gcode_samples, A.gcode_features, B.gcode_features, A.gcode, B.gcode, A.pc.gcode, B.pc.gcode)
colnames(data.gcode) <- c('x', 'y')
data.gcode$matrix <- rep(c('A', 'B', 'A', 'B','A', 'B','A', 'B', 'A', 'B'), each = k)
data.gcode$treatment <- rep(c('1_Original', '2_Encoded Features','3_Encoded Samples', '4_Generative Encoding', '5_Procrustes'), each = 2*k)

library(ggplot2)
data.plot <- ggplot(data.gcode, aes(x = x, y = y, color = matrix)) +
  geom_point(aes(shape = matrix), size = 3, alpha = 0.7) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  facet_wrap(~treatment, scales = 'free') 
data.plot

# ggsave("~/Documents/main_files/AskExplain/gcode/figures/1_heart_procrustes.pdf",data.plot,width = 8,height = 5)

