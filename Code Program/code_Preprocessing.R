setwd("/mgpfs/home/redtrian/Data/convertData")
rm(list=ls())

library(tidyverse)
library(data.table)
library(ggplot2)

# Pastikan sudah meng-autentikasi Google Drive
# download file ZIP manual jika perlu, lalu set path-nya
zip_file_path <- "/mgpfs/home/redtrian/Data/quant.zip"
extract_dir <- "/mgpfs/home/redtrian/Data/convertData/path_to_extract_directory"

# Ekstrak ZIP
unzip(zip_file_path, exdir = extract_dir)

tpm_df <- NULL
file_list <- list.files(extract_dir, pattern = "quant.genes.sf$", recursive = TRUE, full.names = TRUE)

for (file_path in file_list) {
  folder_name <- basename(dirname(file_path))
  df <- fread(file_path, select = c(1, 4), header = TRUE)
  colnames(df) <- c("gene", folder_name)

  if (is.null(tpm_df)) {
    tpm_df <- df
  } else {
    tpm_df <- merge(tpm_df, df, by = "gene", all = TRUE)
  }
}

# Simpan hasil gabungan ke CSV
output_file <- "/mgpfs/home/redtrian/Data/convertData/combined_TPM.csv"
fwrite(tpm_df, output_file)

numerical_cols <- tpm_df[, -1, with = FALSE]  # hilangkan kolom gene
tpm_df$Total <- rowSums(numerical_cols, na.rm = TRUE)

# Simpan ulang
output_file_total <- "/mgpfs/home/redtrian/Data/convertData/combined_TPM_with_total.csv"
fwrite(tpm_df, output_file_total)

tpm_df <- fread(output_file_total)
tpm_df$log2_Total <- log2(tpm_df$Total + 1)

median_tpm <- median(tpm_df$Total, na.rm = TRUE)
median_log2_tpm <- median(tpm_df$log2_Total, na.rm = TRUE)

# Filter berdasarkan median TPM
tpm_df_filtered <- tpm_df[tpm_df$Total >= median_tpm]
fwrite(tpm_df_filtered, "//mgpfs/home/redtrian/Data/convertData/filtered_TPM_for_clustering.csv")

# Buat boxplot
png("/mgpfs/home/redtrian/Data/convertData/tpm_and_log2tpm_boxplots.png", width = 1000, height = 500)
par(mfrow = c(1,2))

boxplot(tpm_df$Total, main = "Boxplot TPM", col = "red", ylab = "TPM")
abline(h = median_tpm, col = "black", lwd = 2)
text(1.2, median_tpm, paste("Median:", round(median_tpm, 2)))

boxplot(tpm_df$log2_Total, main = "Boxplot Log2(TPM+1)", col = "red", ylab = "Log2(TPM+1)")
abline(h = median_log2_tpm, col = "black", lwd = 2)
text(1.2, median_log2_tpm, paste("Median:", round(median_log2_tpm, 2)))

dev.off()

library(scales)

data_filtered <- fread("/mgpfs/home/redtrian/Data/convertData/filtered_TPM_for_clustering.csv")
zscore_data <- data_filtered[, lapply(.SD, scale), .SDcols = !c("gene")]

# Gabungkan kembali dengan nama gen
zscore_data <- cbind(gene = data_filtered$gene, zscore_data)

# Simpan file zscore
fwrite(zscore_data, "/mgpfs/home/redtrian/Data/convertData/zscore_data_with_total.csv")

# Hilangkan kolom Total dan log2_Total
zscore_data_clean <- zscore_data[, !c("Total", "log2_Total"), with = FALSE]
fwrite(zscore_data_clean, "/mgpfs/home/redtrian/Data/convertData/zscore_data.csv")
