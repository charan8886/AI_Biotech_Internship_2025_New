# Install packages (first time only)
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GEOquery")

# Load GEOquery
library(GEOquery)

# Download processed dataset (Series Matrix file)
gse <- getGEO("GSE15852", GSEMatrix = TRUE)

# If multiple platforms, select the first
data <- gse[[1]]

# Expression matrix
expr_data <- exprs(data)

# Metadata (sample info)
pheno_data <- pData(data)

# Save to CSV
write.csv(expr_data, "GSE15852_expression.csv")
write.csv(pheno_data, "GSE15852_metadata.csv")

cat("Processed expression data + metadata saved!\n")



# ----------------------------
# Internship Assignment: Data Retrieval
# Dataset: GSE15852 (Breast Cancer vs Normal Tissue)
# ----------------------------


# ----------------------------
# Step 1: Download dataset from GEO
# ----------------------------
gse_id <- "GSE15852"

# This will download and load the series matrix file
gse <- getGEO(gse_id, GSEMatrix = TRUE)

# Some GEO datasets have multiple platforms, pick the first one
length(gse)   
# check number of platforms
data <- gse[[1]]

# ----------------------------
# Step 2: Explore metadata
# ----------------------------
# Show phenotype data (sample info)
pheno_data <- pData(data)
head(pheno_data)

# Extract titles or characteristics to identify groups
unique(pheno_data$title)

# ----------------------------
# Step 3: Check sample groups
# ----------------------------
# In this dataset, "title" column contains "tumor" or "normal"
table(pheno_data$title)

# Create a simple grouping variable
pheno_data$group <- ifelse(grepl("normal", pheno_data$title, ignore.case = TRUE),
                           "Normal", "Tumor")
table(pheno_data$group)

# ----------------------------
# Step 4: Expression data
# ----------------------------
expr_data <- exprs(data)   # matrix of gene expression values
dim(expr_data)             # rows = genes, columns = samples
expr_data[1:5, 1:5]        # preview first 5x5 section

# ----------------------------
# Step 5: Save for later modules
# ----------------------------
# Save expression matrix
write.csv(expr_data, file = "GSE15852_expression_Charan.csv")

# Save phenotype metadata
write.csv(pheno_data, file = "GSE15852_metadata_Charan.csv")

cat("Dataset download and preprocessing complete. Files saved!\n")









































































