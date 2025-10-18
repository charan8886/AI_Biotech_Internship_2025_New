#############################################################################
# Assignment #4 
#############################################################################

# ---------------------------
# 0) Helper: install & load
# ---------------------------
install_and_load <- function(pkgs)
  {
  for(p in pkgs){
    if(!suppressWarnings(requireNamespace(p, quietly = TRUE))){
      message("Installing package: ", p)
      BiocManager::install(p, ask = FALSE, update = FALSE)
    }
    library(p, character.only = TRUE)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# packages to use
pkgs <- c("GEOquery","limma","AnnotationDbi","AnnotationHub",
          "pheatmap","ggplot2","RColorBrewer","dplyr","Biobase")
install_and_load(pkgs)

# ---------------------------
# 1) Prepare Results folder
# ---------------------------
output_dir <- "Results"
if(!dir.exists(output_dir)) dir.create(output_dir)

# ---------------------------
# 2) Download GEO dataset
# ---------------------------
gse_id <- "GSE15852"
message("Downloading ", gse_id, " (this may take a moment)...")
gse_list <- getGEO(gse_id, GSEMatrix = TRUE)
if(length(gse_list) == 0) stop("No GEO series matrix found for ", gse_id)
eset <- gse_list[[1]]   
# use first platform / ExpressionSet
message("Loaded ExpressionSet with annotation: ", annotation(eset))

# expression matrix and phenotype
expr <- exprs(eset)     # probes x samples
pheno <- pData(eset)    # sample metadata
message("Expression matrix dims (probes x samples): ", paste(dim(expr), collapse = " x "))

# ---------------------------
# 3) Create / verify group labels
# ---------------------------
# Try common columns to create Normal vs Cancer group
if(!"group" %in% colnames(pheno)){
  possible_cols <- c("title","source_name_ch1","characteristics_ch1")
  chosen_col <- NULL
  for(col in possible_cols){
    if(col %in% colnames(pheno)){
      chosen_col <- col; break
    }
  }
  if(is.null(chosen_col)) stop("Cannot find sample descriptor column to create groups. Inspect pData(eset).")
  pheno$group <- ifelse(grepl("normal", pheno[[chosen_col]], ignore.case = TRUE), "Normal", "Cancer")
}
table(pheno$group)

# ---------------------------
# 4) Probe -> Gene mapping (robust attempts)
# ---------------------------
probe_ids <- rownames(expr)
platform_annotation <- annotation(eset)  # e.g., "hgu133plus2" or "GPL570"
message("Platform annotation string: ", platform_annotation)

# Attempt 1: If annotation string looks like GPL570 -> use hgu133plus2.db
anno_pkg_guess <- NA
if(grepl("GPL570", platform_annotation, ignore.case = TRUE) ||
   grepl("hgu133plus2", platform_annotation, ignore.case = TRUE)){
  anno_pkg_guess <- "hgu133plus2.db"
}

# Attempt 2: if annotation string contains "hgu" use that pattern
if(is.na(anno_pkg_guess) && grepl("hgu", platform_annotation, ignore.case = TRUE)){
  anno_pkg_guess <- paste0(tolower(platform_annotation), ".db")
}

mapping_df <- NULL

# Try to install & use guessed annotation package
if(!is.na(anno_pkg_guess)){
  try({
    if(!suppressWarnings(requireNamespace(anno_pkg_guess, quietly = TRUE))){
      BiocManager::install(anno_pkg_guess, ask = FALSE, update = FALSE)
    }
    library(anno_pkg_guess, character.only = TRUE)
    # find keytypes available
    kt <- AnnotationDbi::keytypes(get(anno_pkg_guess))
    # choose common keytype for probes
    keytype_use <- intersect(c("PROBEID","PROBESETID","PROBE"), kt)
    if(length(keytype_use) == 0) keytype_use <- kt[1]
    mapping_df <- AnnotationDbi::select(get(anno_pkg_guess),
                                        keys = probe_ids,
                                        columns = c("SYMBOL","ENTREZID"),
                                        keytype = keytype_use[1])
  }, silent = TRUE)
}

# Fallback 1: Try featureData from the ExpressionSet (some GEO matrices include symbols)
if(is.null(mapping_df) || nrow(mapping_df) == 0){
  fd <- fData(eset)
  if(!is.null(fd)){
    # look for common symbol-like columns
    possible_fd_cols <- c("Gene.symbol","GENE_SYMBOL","Symbol","gene_assignment","Gene Symbol","Gene_symbol")
    found <- intersect(possible_fd_cols, colnames(fd))
    if(length(found) > 0){
      message("Using featureData column: ", found[1], " for mapping.")
      mapping_df <- data.frame(PROBEID = rownames(fd), SYMBOL = fd[[found[1]]], stringsAsFactors = FALSE)
      colnames(mapping_df)[1] <- "ProbeID"
    }
  }
}

# Fallback 2: Use AnnotationHub to try to map via organism genome (slower)
if(is.null(mapping_df) || nrow(mapping_df) == 0){
  message("Attempting AnnotationHub fallback mapping (may take time).")
  try({
    ah <- AnnotationHub::AnnotationHub()
    # try to find an appropriate resource for the platform
    # (This is a fallback and may not always succeed; we keep it optional.)
    # We will not fail here; mapping_df can still be NULL afterwards.
  }, silent = TRUE)
}

if(is.null(mapping_df) || nrow(mapping_df) == 0){
  stop("Probe -> gene symbol mapping failed. Please install the correct annotation package (e.g., hgu133plus2.db) or provide mapping. You can run: BiocManager::install('hgu133plus2.db')")
}

# Standardize mapping frame
# try to detect ProbeID and SYMBOL columns
colnames(mapping_df) <- toupper(colnames(mapping_df))
# lower-case duplicates to find SYMBOL
if("PROBEID" %in% colnames(mapping_df) && "SYMBOL" %in% colnames(mapping_df)){
  mapping_sub <- mapping_df[, c("PROBEID","SYMBOL")]
  colnames(mapping_sub) <- c("ProbeID","SYMBOL")
} else {
  # try other names
  candidates <- colnames(mapping_df)
  pid_col <- candidates[grep("PROBE", candidates, ignore.case = TRUE)][1]
  sym_col <- candidates[grep("SYMBOL|GENE", candidates, ignore.case = TRUE)][1]
  if(is.na(pid_col) || is.na(sym_col)) stop("Unexpected mapping_df structure; please inspect mapping_df.")
  mapping_sub <- mapping_df[, c(pid_col, sym_col)]
  colnames(mapping_sub) <- c("ProbeID","SYMBOL")
}

# Remove empty / NA symbols
mapping_sub$SYMBOL[mapping_sub$SYMBOL == ""] <- NA
mapping_sub <- mapping_sub[!is.na(mapping_sub$SYMBOL), , drop = FALSE]
mapping_sub <- unique(mapping_sub)

message("Mapped probes -> genes (mapped rows): ", nrow(mapping_sub), " out of ", length(probe_ids))

# ---------------------------
# 5) Collapse duplicate probes for same gene
#    Strategy: keep the probe with highest mean expression across samples
# ---------------------------
# compute probe means
probe_means <- rowMeans(expr, na.rm = TRUE)
probe_mean_df <- data.frame(ProbeID = names(probe_means), meanExpr = probe_means, stringsAsFactors = FALSE)

map_mean <- merge(mapping_sub, probe_mean_df, by = "ProbeID", all.x = TRUE)
# remove those with NA SYMBOL or missing mean
map_mean <- map_mean[!is.na(map_mean$SYMBOL) & !is.na(map_mean$meanExpr), ]

# choose probe with maximum mean per SYMBOL
map_mean <- as.data.frame(map_mean)
chosen <- map_mean %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice_max(order_by = meanExpr, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

chosen_probes <- chosen$ProbeID
expr_collapsed <- expr[chosen_probes, , drop = FALSE]
rownames(expr_collapsed) <- chosen$SYMBOL

message("Probes before mapping: ", nrow(expr), "; genes after collapsing: ", nrow(expr_collapsed))

# ---------------------------
# 6) LIMMA differential expression: Cancer vs Normal
# ---------------------------
# ensure pheno rows in same order as columns of expression
if(!all(colnames(expr_collapsed) == rownames(pheno))){
  # reorder pheno to match
  pheno <- pheno[colnames(expr_collapsed), , drop = FALSE]
}

group_factor <- factor(pheno$group)
design <- model.matrix(~0 + group_factor)
colnames(design) <- gsub("group_factor", "", colnames(design))
message("Design matrix columns: ", paste(colnames(design), collapse = ", "))

# check Cancer and Normal present
if(!all(c("Cancer","Normal") %in% levels(group_factor))){
  message("Group levels found: ", paste(levels(group_factor), collapse = ", "))
  stop("Required groups 'Cancer' and 'Normal' not both present. Edit pheno$group accordingly.")
}

contrast.matrix <- makeContrasts(Cancer_vs_Normal = Cancer - Normal, levels = design)

fit <- lmFit(expr_collapsed, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

topTab_all <- topTable(fit2, coef = "Cancer_vs_Normal", number = Inf, sort.by = "P")
deg_all <- data.frame(Gene = rownames(topTab_all),
                      logFC = topTab_all$logFC,
                      AveExpr = topTab_all$AveExpr,
                      P.Value = topTab_all$P.Value,
                      adj.P.Val = topTab_all$adj.P.Val,
                      stringsAsFactors = FALSE)

# thresholds
logfc_cutoff <- 1
padj_cutoff <- 0.05

deg_up <- deg_all[deg_all$logFC > logfc_cutoff & deg_all$adj.P.Val < padj_cutoff, ]
deg_down <- deg_all[deg_all$logFC < -logfc_cutoff & deg_all$adj.P.Val < padj_cutoff, ]

message("DEG counts -> Up: ", nrow(deg_up), " Down: ", nrow(deg_down))

# Save CSVs
write.csv(deg_all, file = file.path(output_dir, "DEG_all_genes_GSE15852.csv"), row.names = FALSE)
write.csv(deg_up, file = file.path(output_dir, "DEG_upregulated_GSE15852.csv"), row.names = FALSE)
write.csv(deg_down, file = file.path(output_dir, "DEG_downregulated_GSE15852.csv"), row.names = FALSE)

# ---------------------------
# 7) Volcano plot (ggplot2) and save PNG
# ---------------------------
volcano_df <- deg_all
volcano_df$Category <- "Not_Significant"
volcano_df$Category[volcano_df$logFC > logfc_cutoff & volcano_df$adj.P.Val < padj_cutoff] <- "Upregulated"
volcano_df$Category[volcano_df$logFC < -logfc_cutoff & volcano_df$adj.P.Val < padj_cutoff] <- "Downregulated"
volcano_df$negLog10P <- -log10(volcano_df$P.Value + 1e-300)

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = negLog10P, color = Category)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not_Significant" = "grey")) +
  theme_minimal() +
  labs(title = paste0("Volcano: Cancer vs Normal (", gse_id, ")"),
       x = "log2 Fold Change",
       y = "-log10(P-value)")

png(file.path(output_dir, "Volcano_GSE15852.png"), width = 1200, height = 900, res = 150)
print(p_volcano)
dev.off()
message("Saved volcano to ", file.path(output_dir, "Volcano_GSE15852.png"))

# ---------------------------
# 8) Heatmap of top 25 DEGs
# ---------------------------
top25 <- head(deg_all[order(deg_all$adj.P.Val), ], 25)
genes_top25 <- top25$Gene
mat_top25 <- expr_collapsed[genes_top25, , drop = FALSE]
# z-score row scaling
mat_top25_z <- t(scale(t(mat_top25)))

sample_annotation <- data.frame(Group = pheno$group)
rownames(sample_annotation) <- rownames(pheno)

png(file.path(output_dir, "Heatmap_Top25_DEGs_GSE15852.png"), width = 1600, height = 1200, res = 150)
pheatmap::pheatmap(mat_top25_z, annotation_col = sample_annotation, show_rownames = TRUE,
                   cluster_cols = TRUE, scale = "row", main = "Top 25 DEGs")
dev.off()
message("Saved heatmap to ", file.path(output_dir, "Heatmap_Top25_DEGs_GSE15852.png"))

# ---------------------------
# 9) Summary text (4-5 lines) saved
# ---------------------------
summary_lines <- c(
  paste0("Dataset: ", gse_id, " (Cancer vs Normal)."),
  "Probe-to-gene mapping: arrays can have multiple probes targeting same gene (different exons/isoforms or redundant probes).",
  "Duplicate handling: for genes with multiple probes we retained the probe with highest mean expression across all samples (choosing the most informative probe).",
  paste0("Contrast: Cancer_vs_Normal. DEGs found with |log2FC| >", logfc_cutoff, " and adj.p <", padj_cutoff, 
         " -> Up: ", nrow(deg_up), ", Down: ", nrow(deg_down), ".")
)
writeLines(summary_lines, con = file.path(output_dir, "DEG_result_summary_GSE15852.txt"))
message("Saved summary to ", file.path(output_dir, "DEG_result_summary_GSE15852.txt"))

# ---------------------------
# 10) List saved files
# ---------------------------
message("All results saved in '", output_dir, "':")
print(list.files(output_dir, full.names = TRUE))

