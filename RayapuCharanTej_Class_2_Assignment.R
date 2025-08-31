install.packages("downloader")
library(downloader)

# 1. BMI_data_1
url_BMI_data_1 <- "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/refs/heads/main/BMI_data_1.csv"


# 2. BMI_data_2
url_BMI_data_2 <- "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/refs/heads/main/BMI_data_2.csv"


# 3. DEGs_data_1 
url_DEGs_data_1 <- "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/refs/heads/main/DEGs_Data_1.csv"


# 4. DEGs_data_2
url_DEGs_data_2 <- "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/refs/heads/main/DEGs_Data_2.csv"


# then create two vwctors, one for URLs and one for file names
file_urls <- c(url_BMI_data_1, url_BMI_data_2, url_DEGs_data_1, url_DEGs_data_2)

file_names <- c("BMI_data_1.csv", "BMI_data_2.csv","DEGs_data_1.csv", "DEGs_data_2.csv")


#now loop through them

for (i in seq_along(file_urls)) {
  download(file_urls[i], file_names[i])
  
}


# download R script using same downloaad function
script_url <- "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/refs/heads/main/Module%20I-Basic_R_Functions-Class_2.R"

script_name <- "Rayapu_charan_tej_class2.R"

download(url = script_url, destfile = script_name)


# open R script (ctrl + o)



# ---------------------------
# AI and Omics Internship - Assignment 2
# ---------------------------

# Step 1: Define classify_gene function
classify_gene <- function(logFC, padj){
  if(is.na(padj)) padj <- 1   # replace NA with 1
  if(logFC > 1 & padj < 0.05){
    return("Upregulated")
  } else if(logFC < -1 & padj < 0.05){
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

# Step 2: Setup input files (GitHub links) and output folder
input_files <- c(
  "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/main/DEGs_Data_1.csv",
  "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/main/DEGs_Data_2.csv"
)

file_names <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

output_dir <- "Results"
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

# Step 3: Process each dataset in loop
result_list <- list()

for(i in seq_along(input_files)){
  cat("\nProcessing:", file_names[i], "\n")
  
  # import dataset directly from GitHub
  data <- read.csv(input_files[i], header = TRUE)
  
  # replace missing padj values with 1
  data$padj[is.na(data$padj)] <- 1
  
  # add 'status' column
  data$status <- mapply(classify_gene, data$logFC, data$padj)
  
  # save processed dataset
  output_file_path <- file.path(output_dir, paste0("Processed_", file_names[i]))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Saved to:", output_file_path, "\n")
  
  # store in results list
  result_list[[file_names[i]]] <- data
  
  # summary counts
  cat("Summary counts:\n")
  print(table(data$status))
}

# Step 4: Save workspace
save.image("Rayapu_Charan_Tej_Class_2_Assignment.RData")
