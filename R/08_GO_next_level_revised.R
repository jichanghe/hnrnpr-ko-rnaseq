# ==========================================================
# NEXT-LEVEL GO / GSEA (clusterProfiler) for:
#   H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/HNRNPR_EMBO.xlsx
#
# Includes:
#   🔥 ORA Up vs Down separated (BP/CC/MF)
#   🎨 EMBO-style dotplots (clean theme, wrapped terms, consistent sizing)
#   📊 Combined BP/CC/MF multi-panel figure (single PDF)
#   🏝 GSEA island / network plot (emapplot) - VERSION-SAFE (no cex_category)
#
# Outputs:
#   H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/GO_nextlevel/
#
# Requirements (install if needed):
#   install.packages(c("readxl","dplyr","tibble","stringr","ggplot2","patchwork"))
#   install.packages("BiocManager")
#   BiocManager::install(c("clusterProfiler","enrichplot","org.Hs.eg.db","aplot"))
# ==========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(stringr)
  
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  
  library(ggplot2)
  library(patchwork)
})

# -------------------- Paths --------------------
infile <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/HNRNPR_EMBO.xlsx"
outdir <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/GO_nextlevel_right"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Settings --------------------
padj_cutoff <- 0.05
lfc_cutoff_for_ORA <- 0.0     # set to 1.0 if you want stricter up/down lists
show_n_each_dir <- 10         # top N terms for up and for down in dotplots
wrap_width <- 45              # wrap long GO names

# -------------------- EMBO-ish theme helpers --------------------
theme_embo <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

wrap_terms <- function(p, width = wrap_width) {
  p + scale_y_discrete(labels = function(x) str_wrap(x, width = width))
}

# -------------------- Read & clean --------------------
df <- read_excel(infile)
df <- as.data.frame(df)

req <- c("gene", "log2FoldChange", "padj")
missing_req <- setdiff(req, colnames(df))
if (length(missing_req) > 0) stop("Missing required columns in Excel: ", paste(missing_req, collapse = ", "))

df2 <- df %>%
  transmute(
    gene  = as.character(gene),
    log2FC = as.numeric(log2FoldChange),
    padj  = as.numeric(padj)
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(log2FC)) %>%
  mutate(gene = str_replace(gene, "\\.\\d+$", "")) %>%   # safety for "GENE.1"
  distinct(gene, .keep_all = TRUE)

# -------------------- Map SYMBOL -> ENTREZ --------------------
conv <- bitr(df2$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
dfm <- df2 %>%
  inner_join(conv, by = c("gene" = "SYMBOL")) %>%
  distinct(ENTREZID, .keep_all = TRUE)

message("Mapped genes (ENTREZ): ", nrow(dfm))
if (nrow(dfm) < 200) {
  message("NOTE: If mapping is unexpectedly low, your gene IDs may be ENSEMBL. ",
          "Change fromType='SYMBOL' to fromType='ENSEMBL'.")
}

universe_entrez <- dfm$ENTREZID

# -------------------- DEG sets for ORA (Up/Down) --------------------
sig_up <- dfm %>%
  filter(!is.na(padj), padj < padj_cutoff, log2FC >  lfc_cutoff_for_ORA) %>%
  pull(ENTREZID)

sig_down <- dfm %>%
  filter(!is.na(padj), padj < padj_cutoff, log2FC < -lfc_cutoff_for_ORA) %>%
  pull(ENTREZID)

message("ORA gene counts: up=", length(sig_up), " down=", length(sig_down))

run_enrichGO <- function(gene_ids, ont) {
  if (length(gene_ids) < 10) return(NULL)
  enrichGO(
    gene          = gene_ids,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont,            # "BP","CC","MF"
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
}

# -------------------- ORA Up/Down for BP/CC/MF --------------------
ora <- list(
  UP = list(
    BP = run_enrichGO(sig_up,   "BP"),
    CC = run_enrichGO(sig_up,   "CC"),
    MF = run_enrichGO(sig_up,   "MF")
  ),
  DOWN = list(
    BP = run_enrichGO(sig_down, "BP"),
    CC = run_enrichGO(sig_down, "CC"),
    MF = run_enrichGO(sig_down, "MF")
  )
)

save_enrich <- function(x, prefix) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(invisible(FALSE))
  write.csv(as.data.frame(x), file.path(outdir, paste0(prefix, ".csv")), row.names = FALSE)
  invisible(TRUE)
}

save_enrich(ora$UP$BP,   "ORA_UP_GO_BP")
save_enrich(ora$UP$CC,   "ORA_UP_GO_CC")
save_enrich(ora$UP$MF,   "ORA_UP_GO_MF")
save_enrich(ora$DOWN$BP, "ORA_DOWN_GO_BP")
save_enrich(ora$DOWN$CC, "ORA_DOWN_GO_CC")
save_enrich(ora$DOWN$MF, "ORA_DOWN_GO_MF")

# -------------------- ORA EMBO-style dotplots --------------------
ora_dot <- function(enrich_obj, title, n = show_n_each_dir) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste0(title, "\n(no significant terms)")))
  }
  p <- dotplot(enrich_obj, showCategory = n) +
    ggtitle(title) +
    theme_embo(11)
  wrap_terms(p, wrap_width)
}

# ORA UP 3-panel
p_ora_up <- (
  ora_dot(ora$UP$BP, "ORA UP: GO Biological Process") |
    ora_dot(ora$UP$CC, "ORA UP: GO Cellular Component") |
    ora_dot(ora$UP$MF, "ORA UP: GO Molecular Function")
) + plot_layout(guides = "collect")

ggsave(file.path(outdir, "ORA_UP_GO_BPCMF_3panel.pdf"),
       p_ora_up, width = 14.5, height = 5.2, device = cairo_pdf)

# ORA DOWN 3-panel
p_ora_down <- (
  ora_dot(ora$DOWN$BP, "ORA DOWN: GO Biological Process") |
    ora_dot(ora$DOWN$CC, "ORA DOWN: GO Cellular Component") |
    ora_dot(ora$DOWN$MF, "ORA DOWN: GO Molecular Function")
) + plot_layout(guides = "collect")

ggsave(file.path(outdir, "ORA_DOWN_GO_BPCMF_3panel.pdf"),
       p_ora_down, width = 14.5, height = 5.2, device = cairo_pdf)

# -------------------- GSEA (BP/CC/MF) --------------------
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

gsea <- list(
  BP = run_gseGO("BP"),
  CC = run_gseGO("CC"),
  MF = run_gseGO("MF")
)

if (!is.null(gsea$BP)) write.csv(as.data.frame(gsea$BP), file.path(outdir, "GSEA_GO_BP.csv"), row.names = FALSE)
if (!is.null(gsea$CC)) write.csv(as.data.frame(gsea$CC), file.path(outdir, "GSEA_GO_CC.csv"), row.names = FALSE)
if (!is.null(gsea$MF)) write.csv(as.data.frame(gsea$MF), file.path(outdir, "GSEA_GO_MF.csv"), row.names = FALSE)

# -------------------- GSEA EMBO-style dotplots (Top Up & Down by NES) --------------------
subset_gsea_top <- function(gsea_obj, n_each = show_n_each_dir) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) return(NULL)
  tab <- as.data.frame(gsea_obj)
  
  up_ids <- tab %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice_head(n = n_each) %>% pull(ID)
  dn_ids <- tab %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice_head(n = n_each) %>% pull(ID)
  keep <- unique(c(up_ids, dn_ids))
  
  sub <- gsea_obj
  sub@result <- gsea_obj@result[gsea_obj@result$ID %in% keep, , drop = FALSE]
  sub
}

gsea_dot <- function(gsea_obj, title) {
  sub <- subset_gsea_top(gsea_obj, show_n_each_dir)
  if (is.null(sub) || nrow(as.data.frame(sub)) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste0(title, "\n(no significant terms)")))
  }
  p <- dotplot(sub, showCategory = nrow(sub@result)) +
    ggtitle(title) +
    theme_embo(11)
  wrap_terms(p, wrap_width)
}

p_gsea_3panel <- (
  gsea_dot(gsea$BP, "GSEA GO BP (Top Up & Down)") |
    gsea_dot(gsea$CC, "GSEA GO CC (Top Up & Down)") |
    gsea_dot(gsea$MF, "GSEA GO MF (Top Up & Down)")
) + plot_layout(guides = "collect")

ggsave(file.path(outdir, "GSEA_GO_BPCMF_3panel_EMBO.pdf"),
       p_gsea_3panel, width = 14.5, height = 6.2, device = cairo_pdf)

# -------------------- 🏝 GSEA “island / network” plot (emapplot) --------------------
# VERSION-SAFE: no cex_category / pie args.
make_emap <- function(gsea_obj, label, n_terms_total = 40, layout = "kk") {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) return(NULL)
  
  tab <- as.data.frame(gsea_obj) %>% arrange(p.adjust)
  keep_ids <- head(tab$ID, n_terms_total)
  
  sub <- gsea_obj
  sub@result <- gsea_obj@result[gsea_obj@result$ID %in% keep_ids, , drop = FALSE]
  if (nrow(sub@result) < 10) return(NULL)
  
  sub_sim <- pairwise_termsim(sub)
  
  p <- tryCatch({
    emapplot(sub_sim, layout = layout, showCategory = nrow(sub_sim@result)) +
      ggtitle(label) +
      theme_embo(11) +
      theme(text = element_text(size = 9))
  }, error = function(e) {
    message("emapplot failed for ", label, ": ", conditionMessage(e))
    NULL
  })
  
  p
}

p_emap_bp <- make_emap(gsea$BP, "GSEA GO BP – network / islands", n_terms_total = 50, layout = "kk")
if (!is.null(p_emap_bp)) {
  ggsave(file.path(outdir, "GSEA_GO_BP_island_network.pdf"),
         p_emap_bp, width = 10.5, height = 8.5, device = cairo_pdf)
} else {
  message("No BP island/network plot produced (not enough significant terms or emapplot failed).")
}

p_emap_cc <- make_emap(gsea$CC, "GSEA GO CC – network / islands", n_terms_total = 40, layout = "kk")
if (!is.null(p_emap_cc)) {
  ggsave(file.path(outdir, "GSEA_GO_CC_island_network.pdf"),
         p_emap_cc, width = 10.5, height = 8.5, device = cairo_pdf)
}

p_emap_mf <- make_emap(gsea$MF, "GSEA GO MF – network / islands", n_terms_total = 40, layout = "kk")
if (!is.null(p_emap_mf)) {
  ggsave(file.path(outdir, "GSEA_GO_MF_island_network.pdf"),
         p_emap_mf, width = 10.5, height = 8.5, device = cairo_pdf)
}

message("\nDONE ✅ Outputs saved in:\n  ", outdir, "\n",
        "Key PDFs:\n",
        "  - ORA_UP_GO_BPCMF_3panel.pdf\n",
        "  - ORA_DOWN_GO_BPCMF_3panel.pdf\n",
        "  - GSEA_GO_BPCMF_3panel_EMBO.pdf\n",
        "  - GSEA_GO_BP_island_network.pdf (and CC/MF if enough terms)\n")
