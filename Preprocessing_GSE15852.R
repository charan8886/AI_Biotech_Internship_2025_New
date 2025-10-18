# AI & Omics Internship Assignment 4
# Preprocessing & Normalization of Microarray Data in R
# Dataset: GSE15852 (Breast Tumor vs Normal)

# ---- Load Required Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "affy", "arrayQualityMetrics", "limma"), ask = FALSE)
install.packages("pheatmap")

library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(limma)
library(pheatmap)

# ---- Download Dataset ----
gse <- getGEO("GSE15852", GSEMatrix = TRUE)
gse <- gse[[1]]

# Extract phenotype (sample) information
pdata <- pData(gse)
head(pdata)

# ---- Quality Control: Before Normalization ----
expr_raw <- exprs(gse)

# Boxplot before normalization
png("Boxplot_Before_Normalization.png", width = 1000, height = 600)
boxplot(expr_raw, main = "Boxplot - Before Normalization", las = 2)
dev.off()

# ---- Normalize the data ----
eset_norm <- normalizeBetweenArrays(expr_raw, method = "quantile")

# Boxplot after normalization
png("Boxplot_After_Normalization.png", width = 1000, height = 600)
boxplot(eset_norm, main = "Boxplot - After Normalization", las = 2)
dev.off()

# ---- Quality Control Report ----
# we can generate a QC report folder in your working directory
arrayQualityMetrics(expressionset = gse, outdir = "QC_Report_GSE15852", force = TRUE)

# ---- Filter low-intensity probes ----
cutoff <- 5
filter <- rowMeans(eset_norm) > cutoff
filtered_data <- eset_norm[filter, ]
cat("Number of probes remaining after filtering:", nrow(filtered_data), "\n")

# ---- PCA Plot ----
pca <- prcomp(t(filtered_data), scale. = TRUE)
group <- as.factor(pdata$source_name_ch1)

png("PCA_After_Normalization.png", width = 1000, height = 600)
plot(pca$x[, 1:2], col = group, pch = 19,
     xlab = "PC1", ylab = "PC2", main = "PCA - After Normalization")
legend("topright", legend = levels(group), col = 1:length(levels(group)), pch = 19)
dev.off()

# ---- Relabel phenotype groups (Normal vs Tumor) ----
pdata$condition <- ifelse(grepl("normal", pdata$source_name_ch1, ignore.case = TRUE),
                          "Normal", "Tumor")

table(pdata$condition)

# ---- Save outputs ----
write.csv(filtered_data, "GSE15852_Normalized_Filtered.csv", row.names = TRUE)
write.csv(pdata, "GSE15852_PhenotypeInfo.csv", row.names = TRUE)




