# ===== PERFECT Heatmap: Top 30 Up + Top 30 Down (RN7SK KO vs WT) =====
# Input:  H:/7SK/7SK/workdir1/results/expression/TPM_EMBO_7SK.xlsx
# Columns: gene, GeneID, RN7SK_WT1-4, RN7SK_KO2-4

library(readxl)
library(dplyr)
library(pheatmap)

# --------- input / output ---------
file_path <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/TPM_EMBO_HNRNPR.xlsx"
out_pdf   <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/TPM_EMBO_HNRNPR_heatmap_top30_UpDown.pdf"
out_png   <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/TPM_EMBO_HNRNPR_heatmap_top30_UpDown.png"

top_n <- 30  # top 30 up + top 30 down

# --------- load data ---------
df <- read_excel(file_path, sheet = 1)
names(df) <- trimws(names(df))  # remove hidden spaces in headers

# --------- sample columns (based on your header) ---------
wt_cols <- c("hnRNPR_WT1","hnRNPR_WT2","hnRNPR_WT3","hnRNPR_WT4")
ko_cols <- c("hnRNPR_KO1","hnRNPR_KO2","hnRNPR_KO3","hnRNPR_KO4")
sample_cols <- c(wt_cols, ko_cols)

missing_cols <- setdiff(sample_cols, names(df))
if (length(missing_cols) > 0) {
  stop("❌ Missing sample columns in Excel:\n  ",
       paste(missing_cols, collapse = ", "),
       "\nCheck with: colnames(read_excel(file_path))")
}

# --------- build row labels: gene (fallback GeneID), unique ---------
gene_sym <- trimws(as.character(df$gene))
gene_id  <- trimws(as.character(df$GeneID))
gene_sym[is.na(gene_sym)] <- ""
gene_id[is.na(gene_id)] <- ""

row_label <- ifelse(nzchar(gene_sym), gene_sym, gene_id)
row_label <- ifelse(nzchar(row_label), row_label, paste0("row_", seq_len(nrow(df))))
row_label <- make.unique(as.character(row_label), sep = "_")

# --------- expression matrix (NO select() masking issues) ---------
expr_df <- as.data.frame(dplyr::select(df, all_of(sample_cols)))
rownames(expr_df) <- row_label

expr <- as.matrix(expr_df)
mode(expr) <- "numeric"

# --------- log2(TPM + 1) ---------
expr_log <- log2(expr + 1)

# --------- KO vs WT mean difference (KO - WT) ---------
WT_mean <- rowMeans(expr_log[, wt_cols, drop = FALSE], na.rm = TRUE)
KO_mean <- rowMeans(expr_log[, ko_cols, drop = FALSE], na.rm = TRUE)
dExpr   <- KO_mean - WT_mean

# --------- select Top Up + Top Down ---------
n_use <- min(top_n, length(dExpr))
top_up <- order(dExpr, decreasing = TRUE)[seq_len(n_use)]
top_dn <- order(dExpr, decreasing = FALSE)[seq_len(n_use)]

keep_idx <- unique(c(top_up, top_dn))
expr_sel <- expr_log[keep_idx, , drop = FALSE]

# --------- row Z-score ---------
expr_z <- t(scale(t(expr_sel)))
expr_z[is.na(expr_z)] <- 0

# --------- column annotation ---------
annotation_col <- data.frame(
  Group = factor(c(rep("WT", length(wt_cols)), rep("KO", length(ko_cols))),
                 levels = c("WT", "KO"))
)
rownames(annotation_col) <- colnames(expr_z)

ann_colors <- list(
  Group = c(WT = "#0072B2", KO = "#D55E00")  # color-blind safe palette
)

# --------- heatmap function ---------
draw_heatmap <- function() {
  pheatmap(
    expr_z,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,   # keep WT -> KO order
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 7,
    fontsize_col = 12,
    border_color = NA,
    main = "Top 30 Up + Top 30 Down (hnRNPR KO vs WT)"
  )
}

# --------- save PDF ---------
pdf(out_pdf, width = 7, height = 9)
draw_heatmap()
dev.off()

# --------- save PNG ---------
png(out_png, width = 2400, height = 3000, res = 300)
draw_heatmap()
dev.off()

message("✅ Heatmap saved to:\n  ", out_pdf, "\n  ", out_png)
