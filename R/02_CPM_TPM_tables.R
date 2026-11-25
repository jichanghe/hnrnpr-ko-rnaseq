# ===== Heatmap: Top 30 Up + Top 30 Down (Ref vs Alt) =====
library(data.table)
library(dplyr)
library(pheatmap)

# --------- input / output ---------
file_path <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/CPM.txt"
out_pdf   <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/CPM_expression_heatmap_top30_UpDown.pdf"
out_png   <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/CPM_expression_heatmap_top30_UpDown.png"

# --------- load data ---------
df <- fread(file_path, data.table = FALSE, sep = "\t", header = TRUE, check.names = FALSE)

# --------- build unique rownames (gene_name preferred) ---------
gene_name <- trimws(as.character(df$gene_name))
gene_id   <- trimws(as.character(df$gene_id))
gene_name[is.na(gene_name)] <- ""

row_id <- ifelse(nzchar(gene_name), gene_name, gene_id)
rownames(df) <- make.unique(row_id, sep = "_")

# --------- sample columns ---------
sample_cols <- c("Ref_1","Ref_2","Ref_3","Ref_4",
                 "Alt_1","Alt_2","Alt_3","Alt_4")

expr <- as.matrix(df[, sample_cols])
mode(expr) <- "numeric"

# --------- log2 transform ---------
expr_log <- log2(expr + 1)

# --------- compute mean Ref / Alt expression ---------
df_means <- data.frame(
  Ref_mean = rowMeans(expr_log[, grep("^Ref", colnames(expr_log))]),
  Alt_mean = rowMeans(expr_log[, grep("^Alt", colnames(expr_log))]),
  dExpr    = rowMeans(expr_log[, grep("^Alt", colnames(expr_log))]) -
    rowMeans(expr_log[, grep("^Ref", colnames(expr_log))])
)

# --------- select Top 30 Up + Top 30 Down ---------
top_up <- order(df_means$dExpr, decreasing = TRUE)[1:30]
top_dn <- order(df_means$dExpr, decreasing = FALSE)[1:30]

keep_idx <- unique(c(top_up, top_dn))
expr_sel <- expr_log[keep_idx, , drop = FALSE]

# --------- row Z-score ---------
expr_z <- t(scale(t(expr_sel)))
expr_z[is.na(expr_z)] <- 0

# --------- column annotation ---------
annotation_col <- data.frame(
  Group = factor(
    ifelse(grepl("^Ref", colnames(expr_z)), "Ref", "Alt"),
    levels = c("Ref", "Alt")
  )
)
rownames(annotation_col) <- colnames(expr_z)

ann_colors <- list(
  Group = c(
    Ref = "#0072B2",   # blue (control)
    Alt = "#D55E00"    # vermillion (perturbation)
  )
)

# --------- draw heatmap ---------
draw_heatmap <- function() {
  pheatmap(
    expr_z,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,   # keep Ref → Alt order (recommended)
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 12,
    border_color = NA,
    main = "Top 30 Up + Top 30 Down Expression (Ref vs Alt)"
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


















# ===== GO enrichment: Top 30 Up vs Down =====
library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# --------- input / output ---------
file_path <- "H:/7SK/7SK/workdir1/results/expression/counts.txt"
out_dir   <- "H:/7SK/7SK/workdir1/results/GO"

dir.create(out_dir, showWarnings = FALSE)

# --------- load expression data ---------
df <- fread(file_path, data.table = FALSE)

# --------- gene names ---------
gene_name <- trimws(as.character(df$gene_name))
gene_id   <- trimws(as.character(df$gene_id))
gene_name[is.na(gene_name)] <- ""

genes <- ifelse(nzchar(gene_name), gene_name, gene_id)

# --------- expression matrix ---------
expr <- as.matrix(df[, c("Ref_1","Ref_2","Ref_3","Ref_4",
                         "Alt_1","Alt_2","Alt_3","Alt_4")])
mode(expr) <- "numeric"
expr_log <- log2(expr + 1)

# --------- compute Ref vs Alt change ---------
Ref_mean <- rowMeans(expr_log[, grep("^Ref", colnames(expr_log))])
Alt_mean <- rowMeans(expr_log[, grep("^Alt", colnames(expr_log))])
dExpr    <- Alt_mean - Ref_mean

# --------- select Top 30 Up & Down ---------
top_up_genes <- genes[order(dExpr, decreasing = TRUE)][1:30]
top_dn_genes <- genes[order(dExpr, decreasing = FALSE)][1:30]

# --------- remove duplicates ---------
top_up_genes <- unique(top_up_genes)
top_dn_genes <- unique(top_dn_genes)

# --------- SYMBOL -> ENTREZ ---------
up_entrez <- bitr(
  top_up_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)$ENTREZID

dn_entrez <- bitr(
  top_dn_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)$ENTREZID

# --------- GO Biological Process ---------
ego_up <- enrichGO(
  gene          = up_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_dn <- enrichGO(
  gene          = dn_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# --------- save result tables ---------
write.csv(as.data.frame(ego_up),
          file.path(out_dir, "GO_BP_Top30_Up.csv"),
          row.names = FALSE)

write.csv(as.data.frame(ego_dn),
          file.path(out_dir, "GO_BP_Top30_Down.csv"),
          row.names = FALSE)

# --------- dotplots ---------
pdf(file.path(out_dir, "GO_BP_Top30_Up_dotplot.pdf"), width = 7, height = 5)
dotplot(ego_up, showCategory = 15, title = "GO BP – Top 30 Up genes")
dev.off()

pdf(file.path(out_dir, "GO_BP_Top30_Down_dotplot.pdf"), width = 7, height = 5)
dotplot(ego_dn, showCategory = 15, title = "GO BP – Top 30 Down genes")
dev.off()

message("✅ GO analysis complete. Results saved in: ", out_dir)




