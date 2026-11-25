# ============================================================
# clusterProfiler GO (BP/CC/MF) from DE table Excel:
#   H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/HNRNPR_EMBO.xlsx
#
# Uses columns:
#   gene (SYMBOL preferred), log2FoldChange, padj (or pvalue)
#
# Outputs:
#   - ORA GO dotplots (BP/CC/MF) + CSV tables
#   - GSEA GO dotplots (BP/CC/MF) + CSV tables (recommended)
#
# NOTE:
#   - GO enrichment should use DE results (log2FC, padj), NOT VST/TPM/CPM.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(stringr)
  
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  
  library(ggplot2)
})

# ----------------- INPUT / OUTPUT -----------------
infile <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/HNRNPR_EMBO.xlsx"
outdir <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/GO_clusterProfiler"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sheet_name <- NULL  # set to "Sheet1" etc if needed

# Thresholds for ORA (cutoff-based enrichment)
padj_cutoff <- 0.05
lfc_cutoff  <- 0.0   # set to 1 if you want |log2FC|>1 for ORA lists

# Plot / result settings
show_n_terms <- 20
set.seed(1)

# ----------------- READ EXCEL -----------------
df <- if (is.null(sheet_name)) read_excel(infile) else read_excel(infile, sheet = sheet_name)
df <- as.data.frame(df)

# ----------------- COLUMN DETECTION -----------------
# Prefer gene symbol in column "gene"; otherwise try "genename" or "Gene" or "ID"
gene_col_candidates <- c("gene", "genename", "Gene", "SYMBOL", "symbol", "ID")
gene_col <- gene_col_candidates[gene_col_candidates %in% colnames(df)][1]
if (is.na(gene_col)) stop("Could not find a gene column. Tried: ", paste(gene_col_candidates, collapse = ", "))

# log2FC column
lfc_col_candidates <- c("log2FoldChange", "logFC", "log2FC", "log2_fold_change")
lfc_col <- lfc_col_candidates[lfc_col_candidates %in% colnames(df)][1]
if (is.na(lfc_col)) stop("Could not find log2FC column. Tried: ", paste(lfc_col_candidates, collapse = ", "))

# padj column (fallback to pvalue if padj missing)
padj_col <- if ("padj" %in% colnames(df)) "padj" else if ("p.adjust" %in% colnames(df)) "p.adjust" else NA
pval_col <- if ("pvalue" %in% colnames(df)) "pvalue" else if ("PValue" %in% colnames(df)) "PValue" else NA

if (is.na(padj_col) && is.na(pval_col)) {
  stop("Could not find padj/p.adjust OR pvalue/PValue in the table.")
}

# ----------------- CLEAN TABLE -----------------
df2 <- df %>%
  transmute(
    SYMBOL = as.character(.data[[gene_col]]),
    log2FC = as.numeric(.data[[lfc_col]]),
    padj   = if (!is.na(padj_col)) as.numeric(.data[[padj_col]]) else NA_real_,
    pvalue = if (!is.na(pval_col)) as.numeric(.data[[pval_col]]) else NA_real_
  ) %>%
  filter(!is.na(SYMBOL), SYMBOL != "", !is.na(log2FC)) %>%
  mutate(
    # common cleanup: remove version suffixes if present like "GENE.1"
    SYMBOL = str_replace(SYMBOL, "\\.\\d+$", "")
  ) %>%
  distinct(SYMBOL, .keep_all = TRUE)

# If padj missing, warn & use pvalue in ORA (not ideal)
use_padj <- !all(is.na(df2$padj))
if (!use_padj) message("WARNING: padj not found; ORA will use pvalue cutoff instead (less stringent).")

# ----------------- MAP SYMBOL -> ENTREZ -----------------
# clusterProfiler works best with ENTREZID
conv <- bitr(df2$SYMBOL,
             fromType = "SYMBOL",
             toType   = "ENTREZID",
             OrgDb    = org.Hs.eg.db)

# merge mapping back
dfm <- df2 %>%
  inner_join(conv, by = c("SYMBOL" = "SYMBOL")) %>%
  distinct(ENTREZID, .keep_all = TRUE)

if (nrow(dfm) < 100) {
  message("Mapped ENTREZ genes: ", nrow(dfm),
          " (If this is unexpectedly low, check if your gene IDs are SYMBOL or ENSEMBL.)")
}

# Define universe/background as all tested genes (recommended)
universe_entrez <- dfm$ENTREZID

# ----------------- ORA GENE LISTS (padj cutoff) -----------------
if (use_padj) {
  sig_up   <- dfm %>% filter(!is.na(padj), padj < padj_cutoff, log2FC >  lfc_cutoff) %>% pull(ENTREZID)
  sig_down <- dfm %>% filter(!is.na(padj), padj < padj_cutoff, log2FC < -lfc_cutoff) %>% pull(ENTREZID)
  sig_all  <- dfm %>% filter(!is.na(padj), padj < padj_cutoff, abs(log2FC) > lfc_cutoff) %>% pull(ENTREZID)
} else {
  sig_up   <- dfm %>% filter(!is.na(pvalue), pvalue < padj_cutoff, log2FC >  lfc_cutoff) %>% pull(ENTREZID)
  sig_down <- dfm %>% filter(!is.na(pvalue), pvalue < padj_cutoff, log2FC < -lfc_cutoff) %>% pull(ENTREZID)
  sig_all  <- dfm %>% filter(!is.na(pvalue), pvalue < padj_cutoff, abs(log2FC) > lfc_cutoff) %>% pull(ENTREZID)
}

message("ORA genes: all=", length(sig_all), " up=", length(sig_up), " down=", length(sig_down))

# ----------------- ORA: enrichGO function -----------------
run_enrichGO <- function(gene_ids, ont) {
  if (length(gene_ids) < 10) return(NULL)
  enrichGO(
    gene          = gene_ids,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,            # "BP" "CC" "MF"
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
}

# ORA for ALL sig genes (you can also do UP/DOWN separately if you want)
ego_bp <- run_enrichGO(sig_all, "BP")
ego_cc <- run_enrichGO(sig_all, "CC")
ego_mf <- run_enrichGO(sig_all, "MF")

# Save ORA tables
save_enrich_table <- function(enrich_obj, outfile_prefix) {
  if (is.null(enrich_obj)) return(invisible(FALSE))
  tab <- as.data.frame(enrich_obj)
  write.csv(tab, file.path(outdir, paste0(outfile_prefix, ".csv")), row.names = FALSE)
  invisible(TRUE)
}

save_enrich_table(ego_bp, "ORA_GO_BP")
save_enrich_table(ego_cc, "ORA_GO_CC")
save_enrich_table(ego_mf, "ORA_GO_MF")

# ORA dotplots
plot_dot <- function(enrich_obj, title, outfile) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) return(invisible(FALSE))
  p <- dotplot(enrich_obj, showCategory = show_n_terms) +
    ggtitle(title) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_blank()
    )
  ggsave(file.path(outdir, outfile), p, width = 8.5, height = 5.5, device = cairo_pdf)
  invisible(TRUE)
}

plot_dot(ego_bp, "ORA GO: Biological Process", "ORA_GO_BP_dotplot.pdf")
plot_dot(ego_cc, "ORA GO: Cellular Component", "ORA_GO_CC_dotplot.pdf")
plot_dot(ego_mf, "ORA GO: Molecular Function", "ORA_GO_MF_dotplot.pdf")

# ----------------- GSEA (RECOMMENDED) -----------------
# Rank ALL genes by log2FC (named numeric vector: ENTREZID -> log2FC)
geneList <- dfm$log2FC
names(geneList) <- dfm$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

run_gseGO <- function(ont) {
  gseGO(
    geneList      = geneList,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,
    minGSSize     = 10,
    maxGSSize     = 500,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    verbose       = FALSE
  )
}

gsea_bp <- run_gseGO("BP")
gsea_cc <- run_gseGO("CC")
gsea_mf <- run_gseGO("MF")

# Save GSEA tables
if (!is.null(gsea_bp)) write.csv(as.data.frame(gsea_bp), file.path(outdir, "GSEA_GO_BP.csv"), row.names = FALSE)
if (!is.null(gsea_cc)) write.csv(as.data.frame(gsea_cc), file.path(outdir, "GSEA_GO_CC.csv"), row.names = FALSE)
if (!is.null(gsea_mf)) write.csv(as.data.frame(gsea_mf), file.path(outdir, "GSEA_GO_MF.csv"), row.names = FALSE)

# GSEA dotplots (split by direction using NES)
plot_gsea_dot <- function(gsea_obj, title, outfile) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) return(invisible(FALSE))
  tab <- as.data.frame(gsea_obj)
  tab$.dir <- ifelse(tab$NES >= 0, "Activated (NES>0)", "Suppressed (NES<0)")
  # show top terms from each direction
  tab2 <- bind_rows(
    tab %>% filter(.dir == "Activated (NES>0)") %>% arrange(p.adjust) %>% slice_head(n = show_n_terms),
    tab %>% filter(.dir == "Suppressed (NES<0)") %>% arrange(p.adjust) %>% slice_head(n = show_n_terms)
  )
  # Make dotplot from filtered result by subsetting gsea object
  keep_ids <- tab2$ID
  gsea_sub <- gsea_obj
  gsea_sub@result <- gsea_obj@result[gsea_obj@result$ID %in% keep_ids, , drop = FALSE]
  
  p <- dotplot(gsea_sub, showCategory = nrow(gsea_sub@result), split = ".sign") +
    ggtitle(title) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_blank()
    )
  
  ggsave(file.path(outdir, outfile), p, width = 8.5, height = 7.0, device = cairo_pdf)
  invisible(TRUE)
}

# If you prefer a simpler plot, you can just do dotplot(gsea_bp, showCategory=20)
plot_gsea_dot(gsea_bp, "GSEA GO: Biological Process", "GSEA_GO_BP_dotplot.pdf")
plot_gsea_dot(gsea_cc, "GSEA GO: Cellular Component", "GSEA_GO_CC_dotplot.pdf")
plot_gsea_dot(gsea_mf, "GSEA GO: Molecular Function", "GSEA_GO_MF_dotplot.pdf")

# ----------------- DONE -----------------
message("DONE ✅ Results saved in: ", outdir)

# Helpful note:
# If SYMBOL mapping is poor (many genes fail), your gene IDs may be ENSEMBL.
# In that case replace the mapping step with:
#   bitr(df2$SYMBOL, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
