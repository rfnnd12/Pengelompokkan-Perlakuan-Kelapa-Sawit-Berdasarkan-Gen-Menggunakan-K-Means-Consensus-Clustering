# ================================================================
# Langsung Menghitung Silhouette untuk Hasil Klasterisasi ConsensusClusterPlus
# ================================================================
# Set working directory
setwd("/mgpfs/home/redtrian/Data/Demo_semhas")
rm(list=ls())

library(cluster)
library(factoextra)
library(data.table)

load("/mgpfs/home/redtrian/Data/Demo_semhas/output/ConsCluster_maxK50_rep20.RData")

# Membaca data
data <- data.frame(fread("/mgpfs/home/redtrian/Data/Demo_semhas/zscore_data.csv"))
# Membuat kolom pertama menjadi nama baris
rownames(data) <- data$gene
# Menghilangkan kolom pertama
data <- data[,-1]
data <- as.matrix(data)

silhouette_widths <- numeric()  
maxK <- 50
# Loop untuk menghitung silhouette untuk setiap jumlah klaster k dari hasil ConsensusClusterPlus
for (k in 2:maxK) {

  # Ambil hasil klasterisasi untuk k clusters dari ConsensusClusterPlus
  kmeans_result <- results[[k]]$consensusClass
  
  # Hitung jarak Euclidean antar data
  dist_matrix <- dist(t(data))
  
  # Pastikan dist_matrix adalah matriks yang kompatibel dengan fungsi silhouette
  dist_matrix <- as.dist(dist_matrix)
  
  # Hitung nilai silhouette menggunakan hasil klasterisasi dan matriks jarak
  sil <- silhouette(kmeans_result, dist_matrix)
  
  # Simpan rata-rata nilai silhouette
  silhouette_widths[k] <- mean(sil[, 3])
}

# Simpan plot Silhouette dalam format PNG
outFolder <- "/mgpfs/home/redtrian/Data/Demo_semhas/silhouette/" 
dir.create(outFolder, recursive = TRUE)
silhouette_plot <- paste0(outFolder, "/silhouette_plot_consensus_k", maxK, ".png")  
png(silhouette_plot, width = 800, height = 600)

# Plot nilai Average Silhouette Width
plot(2:maxK, silhouette_widths[2:maxK], type = "b", pch = 19,
     frame = FALSE, xlab = "Number of clusters k", ylab = "Average silhouette width",
     main = "Optimal number of clusters (ConsensusClusterPlus)")

# Menambahkan garis vertikal untuk jumlah klaster optimal
optimal_k <- which.max(silhouette_widths)
abline(v = optimal_k, lty = 2, col = "blue")

# Selesai menyimpan plot
dev.off()

# Menampilkan jumlah klaster optimal dan rata-rata Silhouette Width
cat("Optimal number of clusters (k):", optimal_k, "\n")
cat("Average Silhouette Width for optimal k:", silhouette_widths[optimal_k], "\n")

# Menyimpan hasil rata-rata Silhouette Width untuk setiap k ke file teks
output_file <- "/mgpfs/home/redtrian/Data/Demo_semhas/silhouette/silhouette_widths.txt"
write.table(silhouette_widths[2:maxK], file = output_file, col.names = FALSE, row.names = TRUE, quote = FALSE)

cat("Rata-rata Silhouette Width telah disimpan ke:", output_file, "\n")