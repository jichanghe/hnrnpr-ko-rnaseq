# ============================================================
# VST heatmap from Excel:
#   H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/HNRNPR_EMBO.xlsx
#
# Uses raw counts columns:
#   hnRNPR_KO1-4, hnRNPR_WT1-4
# Runs DESeq2 -> vst(blind=FALSE)
# Exports Illustrator-ready vector PDF heatmap
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(DESeq2)
  library(matrixStats)
  library(pheatmap)
  library(tibble)
  library(dplyr)
  library(stringr)
})

# --------- INPUT / OUTPUT ----------
infile  <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/HNRNPR_EMBO.xlsx"
out_pdf <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/HNRNPR_VST_heatmap.pdf"

# If your Excel has multiple sheets, set it here (NULL = first sheet)
sheet_name <- NULL

# Heatmap selection mode:
#   "padj"     = use padj < padj_cutoff (falls back to top variable if too few)
#   "top_var"  = ignore padj; use top N most variable genes
selection_mode <- "padj"

padj_cutoff <- 0.05
min_genes_if_padj <- 30   # if fewer than this pass padj, fall back to top variable
top_n_var <- 50           # used for "top_var" or fallback

# --------- COLUMNS ----------
count_cols <- c(
  "hnRNPR_KO1","hnRNPR_KO2","hnRNPR_KO3","hnRNPR_KO4",
  "hnRNPR_WT1","hnRNPR_WT2","hnRNPR_WT3","hnRNPR_WT4"
)

# Prefer gene symbol column name "gene"; if missing, try "ID"
gene_col_candidates <- c("gene", "ID")

# --------- READ EXCEL ----------
df <- if (is.null(sheet_name)) {
  read_excel(infile)
} else {
  read_excel(infile, sheet = sheet_name)
}

df <- as.data.frame(df)

# --------- BASIC CHECKS ----------
missing_counts <- setdiff(count_cols, colnames(df))
if (length(missing_counts) > 0) {
  stop("Missing count columns in Excel: ", paste(missing_counts, collapse = ", "))
}

gene_col <- gene_col_candidates[gene_col_candidates %in% colnames(df)][1]
if (is.na(gene_col) || is.null(gene_col)) {
  stop("Could not find a gene ID column. Expected one of: ", paste(gene_col_candidates, collapse = ", "))
}

# Make sure counts are numeric integers
for (cc in count_cols) {
  df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
}

# Remove rows with no gene name
df <- df[!is.na(df[[gene_col]]) & df[[gene_col]] != "", , drop = FALSE]

# Replace NA counts with 0 (rare, but safer)
for (cc in count_cols) df[[cc]][is.na(df[[cc]])] <- 0

# Force non-negative counts
for (cc in count_cols) {
  if (any(df[[cc]] < 0, na.rm = TRUE)) stop("Found negative counts in column: ", cc)
}

# Deduplicate genes:
# keep the row with the largest total counts across samples
df <- df %>%
  mutate(.gene = as.character(.data[[gene_col]]),
         .row_sum = rowSums(across(all_of(count_cols)), na.rm = TRUE)) %>%
  arrange(desc(.row_sum)) %>%
  distinct(.gene, .keep_all = TRUE)

# Build count matrix (genes x samples)
count_mat <- as.matrix(df[, count_cols, drop = FALSE])
rownames(count_mat) <- df$.gene

# Ensure integer mode for DESeq2
storage.mode(count_mat) <- "integer"

# Filter extremely low-count genes (helps speed + stability)
keep <- rowSums(count_mat) >= 10
count_mat <- count_mat[keep, , drop = FALSE]

if (nrow(count_mat) < 10) {
  stop("Too few genes left after filtering. Check your input counts.")
}

# --------- METADATA ----------
coldata <- data.frame(
  condition = factor(c("KO","KO","KO","KO","WT","WT","WT","WT"), levels = c("WT","KO"))
)
rownames(coldata) <- colnames(count_mat)

# --------- DESEQ2 + VST ----------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

# --------- PICK GENES FOR HEATMAP ----------
use_mat <- NULL

if (selection_mode == "padj") {
  # Use padj from your Excel file if present; otherwise compute DESeq2 results
  if ("padj" %in% colnames(df)) {
    padj_map <- df %>% select(.gene, padj)
    padj_vec <- padj_map$padj
    names(padj_vec) <- padj_map$.gene
    
    padj_vec <- padj_vec[rownames(vst_mat)]
    sig <- names(padj_vec)[!is.na(padj_vec) & padj_vec < padj_cutoff]
  } else {
    res <- results(dds)
    sig <- rownames(res)[!is.na(res$padj) & res$padj < padj_cutoff]
  }
  
  if (length(sig) >= min_genes_if_padj) {
    use_mat <- vst_mat[sig, , drop = FALSE]
  } else {
    # fallback to top variable genes
    vars <- rowVars(vst_mat)
    top  <- head(order(vars, decreasing = TRUE), top_n_var)
    use_mat <- vst_mat[top, , drop = FALSE]
    message("Note: too few genes passed padj cutoff; using top variable genes instead (n=", nrow(use_mat), ").")
  }
  
} else if (selection_mode == "top_var") {
  vars <- rowVars(vst_mat)
  top  <- head(order(vars, decreasing = TRUE), top_n_var)
  use_mat <- vst_mat[top, , drop = FALSE]
} else {
  stop("selection_mode must be 'padj' or 'top_var'")
}

# Optional: cap extreme values for nicer color scaling (helps outliers)
cap <- 3
use_scaled <- t(scale(t(use_mat)))                # row Z-score
use_scaled[use_scaled >  cap] <-  cap
use_scaled[use_scaled < -cap] <- -cap
use_scaled[is.na(use_scaled)] <- 0

# --------- ANNOTATION BAR (WT/KO) ----------
ann_col <- data.frame(condition = coldata$condition)
rownames(ann_col) <- rownames(coldata)

# --------- EXPORT HEATMAP PDF (vector) ----------
# pheatmap draws as vector in pdf() (Illustrator-friendly)
pdf(out_pdf, width = 7.2, height = 9.0, useDingbats = FALSE)

pheatmap(
  use_scaled,
  annotation_col = ann_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 10,
  border_color = NA,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete"
)

dev.off()

message("DONE ✅ Heatmap saved to: ", out_pdf)

