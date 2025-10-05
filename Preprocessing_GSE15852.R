
# Preprocessing and Normalization of Microarray Data
# Dataset: GSE15852
# Author: Rayapu Charan Tej

# ============================
# 1 Load required packages
# ============================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("affy", "oligo", "arrayQualityMetrics", "GEOquery", "limma"), ask = FALSE)

library(affy)
library(oligo)
library(arrayQualityMetrics)
library(GEOquery)
library(limma)

# ============================
# 2️ Download and load the dataset
# ============================
# Get the GEO dataset
gse <- getGEO("GSE15852", GSEMatrix = TRUE)
exprSet <- exprs(gse[[1]])        # Expression matrix
pdata <- pData(gse[[1]])          # Phenotype / sample info

# Check data dimensions
dim(exprSet)
head(pdata[,1:5])

# ============================
# 3 Quality control 
# ============================

# Boxplot before normalization
boxplot(exprSet, main = "Boxplot - Before Normalization", las = 2)

# Sample correlation heatmap
library(pheatmap)
cor_matrix <- cor(exprSet)
pheatmap(cor_matrix, main = "Sample Correlation - Before Normalization")

# Identify potential outliers (based on distance)
sampleDists <- dist(t(exprSet))
plot(hclust(sampleDists), main = "Hierarchical Clustering - Before Normalization")

#  Save QC report
arrayQualityMetrics(
  expressionset = gse[[1]],
  outdir = "QC_Before_Normalization",
  force = TRUE,
  do.logtransform = TRUE
)

# ============================
# 4️ Normalization
# ============================

# RMA normalization (Robust Multi-array Average)
eset_norm <- normalizeBetweenArrays(exprSet, method = "quantile")

# Boxplot after normalization
boxplot(eset_norm, main = "Boxplot - After Normalization", las = 2)

# QC again after normalization
arrayQualityMetrics(
  ExpressionSet(eset_norm, phenoData = AnnotatedDataFrame(pdata)),
  outdir = "QC_After_Normalization",
  force = TRUE,
  do.logtransform = FALSE
)

# ============================
# 5️ Filtering low-intensity probes
# ============================
# Remove probes with low expression across samples
threshold <- 5
keep <- rowMeans(eset_norm) > threshold
filtered_data <- eset_norm[keep, ]

cat("Number of probes before filtering:", nrow(eset_norm), "\n")
cat("Number of probes after filtering:", nrow(filtered_data), "\n")

# ============================
# 6️ Re-label phenotype data (Normal vs Tumor)
# ============================

# Inspect sample titles to identify groups
unique(pdata$title)

# Create condition label based on title
pdata$condition <- ifelse(grepl("normal", pdata$title, ignore.case = TRUE),
                          "Normal", "Tumor")

table(pdata$condition)

# ============================
# 7️ Save processed data
# ============================
# Combine filtered expression matrix and phenotype info
write.csv(filtered_data, "GSE15852_Normalized_Filtered.csv")
write.csv(pdata, "GSE15852_PhenotypeInfo.csv")



