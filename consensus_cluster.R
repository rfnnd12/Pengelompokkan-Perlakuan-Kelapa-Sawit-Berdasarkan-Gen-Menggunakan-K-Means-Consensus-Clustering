setwd("/mgpfs/home/redtrian/Data/clus50")
rm(list=ls())

library(data.table)
library(ConsensusClusterPlus)

#membaca data
data <- data.frame(fread("/mgpfs/home/redtrian/Data/zscore_data.csv"))
#membuat kolom pertama menjadi nama baris
rownames(data) <- data$gene
#hilangkan kolom 1
data <- data[,-1]
data <- as.matrix(data)

outFolder <- "/mgpfs/home/redtrian/Data/clus50/output/"
maxK <- 50
rep <- 20

results <- ConsensusClusterPlus(data, maxK=maxK, reps=rep, pItem=0.8, pFeature=1, title ="output",
                               clusterAlg="km", distance="euclidean",
                               seed=12345, plot="pdf", writeTable=TRUE)

save(results, file = paste0(outFolder, "/ConsCluster_maxK", maxK, "_rep", rep, ".RData"))
