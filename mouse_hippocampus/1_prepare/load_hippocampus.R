setwd("")
expression <- read.table("./data/external/single_nucleus_rnaseq_adult_mouse_hippocampus/DATA_MATRIX_LOG_TPM.txt",header=T,row.names = 1)
cell_clusters <- read.table("./data/external/single_nucleus_rnaseq_adult_mouse_hippocampus/CLUSTER_AND_SUBCLUSTER_INDEX.txt",header=T,skip=1)
