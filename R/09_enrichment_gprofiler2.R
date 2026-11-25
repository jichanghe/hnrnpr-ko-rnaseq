# ==========================================================
# EMBO / Nature-ready enrichment figure set:
#   ✅ GSEA GO BP NES-centered (Up on top, Down bottom)
#   ✅ Color by NES (red=Activated, blue=Suppressed)
#   ✅ Numeric NES labels
#   ✅ Clearer GSEA island/network (term similarity + clustering)
#   ✅ Add KEGG + Reactome in same EMBO layout (GSEA + optional ORA)
#
# Input:
#   H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/HNRNPR_EMBO.xlsx
#
# Output folder:
#   H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/ENRICH_EMBO/
#
# Install (if needed):
#   install.packages(c("readxl","dplyr","tibble","stringr","ggplot2","patchwork","ggrepel"))
#   install.packages("BiocManager")
#   BiocManager::install(c("clusterProfiler","enrichplot","org.Hs.eg.db","ReactomePA"))
# ==========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(stringr)
  
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ReactomePA)
  
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

# -------------------- Paths --------------------
infile <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/HNRNPR_EMBO.xlsx"
outdir <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/ENRICH_EMBO"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Settings --------------------
padj_cutoff <- 0.05
wrap_width <- 45

# GSEA selection for dotplots
topN_each_dir <- 10       # top NES>0 and top NES<0
min_terms_for_network <- 15
network_keep_terms <- 60  # how many top terms (by p.adjust) to keep for island/network

# KEGG organism (human)
kegg_org <- "hsa"

# -------------------- EMBO theme helpers --------------------
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

wrap_y <- function(p, width = wrap_width) {
  p + scale_y_discrete(labels = function(x) str_wrap(x, width = width))
}

# Safe save
save_pdf <- function(plot, filename, w, h) {
  ggsave(file.path(outdir, filename), plot, width = w, height = h, device = cairo_pdf)
}

# -------------------- Read & clean --------------------
df <- read_excel(infile) |> as.data.frame()

req <- c("gene", "log2FoldChange", "padj")
miss <- setdiff(req, colnames(df))
if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

df2 <- df |>
  transmute(
    gene  = as.character(gene),
    log2FC = as.numeric(log2FoldChange),
    padj  = as.numeric(padj)
  ) |>
  filter(!is.na(gene), gene != "", !is.na(log2FC)) |>
  mutate(gene = str_replace(gene, "\\.\\d+$", "")) |>
  distinct(gene, .keep_all = TRUE)

# -------------------- Map SYMBOL -> ENTREZ --------------------
conv <- bitr(df2$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
dfm <- df2 |>
  inner_join(conv, by = c("gene" = "SYMBOL")) |>
  distinct(ENTREZID, .keep_all = TRUE)

message("Mapped genes (ENTREZ): ", nrow(dfm))
if (nrow(dfm) < 200) {
  message("NOTE: If mapping is unexpectedly low, your gene IDs may be ENSEMBL. ",
          "Change fromType='SYMBOL' to fromType='ENSEMBL'.")
}

# -------------------- Ranked list for GSEA --------------------
geneList <- dfm$log2FC
names(geneList) <- dfm$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

# -------------------- GSEA runners --------------------
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

run_gseKEGG <- function() {
  gseKEGG(
    geneList      = geneList,
    organism      = kegg_org,
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
}

run_gseReactome <- function() {
  ReactomePA::gsePathway(
    geneList      = geneList,
    organism      = "human",
    minGSSize     = 10,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
}

# -------------------- Compute enrichments --------------------
gsea_go_bp <- run_gseGO("BP")
gsea_go_cc <- run_gseGO("CC")
gsea_go_mf <- run_gseGO("MF")

gsea_kegg <- run_gseKEGG()
gsea_react <- run_gseReactome()

# Save tables
if (!is.null(gsea_go_bp)) write.csv(as.data.frame(gsea_go_bp), file.path(outdir, "GSEA_GO_BP.csv"), row.names = FALSE)
if (!is.null(gsea_go_cc)) write.csv(as.data.frame(gsea_go_cc), file.path(outdir, "GSEA_GO_CC.csv"), row.names = FALSE)
if (!is.null(gsea_go_mf)) write.csv(as.data.frame(gsea_go_mf), file.path(outdir, "GSEA_GO_MF.csv"), row.names = FALSE)
if (!is.null(gsea_kegg))  write.csv(as.data.frame(gsea_kegg),  file.path(outdir, "GSEA_KEGG.csv"),  row.names = FALSE)
if (!is.null(gsea_react)) write.csv(as.data.frame(gsea_react), file.path(outdir, "GSEA_Reactome.csv"), row.names = FALSE)

# ==========================================================
# 1) NES-centered dotplot function (Up top / Down bottom)
#    - y axis is term
#    - x axis is NES
#    - color reflects NES sign (blue->red)
#    - size = setSize
#    - label = numeric NES
# ==========================================================
prep_top_up_down <- function(gsea_obj, n_each = topN_each_dir) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) return(NULL)
  tab <- as.data.frame(gsea_obj)
  
  up <- tab |> filter(NES > 0) |> arrange(p.adjust) |> slice_head(n = n_each)
  dn <- tab |> filter(NES < 0) |> arrange(p.adjust) |> slice_head(n = n_each)
  
  tab2 <- bind_rows(
    up |> mutate(.dir = "Activated (NES>0)"),
    dn |> mutate(.dir = "Suppressed (NES<0)")
  )
  
  # Order so Activated on top, Suppressed bottom, and within each by NES magnitude / significance
  tab2 <- tab2 |>
    mutate(Description = str_wrap(Description, width = wrap_width)) |>
    mutate(.dir = factor(.dir, levels = c("Activated (NES>0)", "Suppressed (NES<0)"))) |>
    arrange(.dir, desc(NES)) |>
    mutate(Description = factor(Description, levels = rev(unique(Description))))
  
  tab2
}

nes_centered_plot <- function(gsea_obj, title) {
  tab2 <- prep_top_up_down(gsea_obj, topN_each_dir)
  if (is.null(tab2) || nrow(tab2) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste0(title, "\n(no significant terms)")))
  }
  
  ggplot(tab2, aes(x = NES, y = Description)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
    geom_point(aes(size = setSize, color = NES)) +
    ggrepel::geom_text_repel(
      aes(label = sprintf("%.2f", NES)),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.2,
      point.padding = 0.15,
      segment.size = 0.2,
      min.segment.length = 0
    ) +
    scale_color_gradient2(
      low = "#2b6cb0", mid = "grey90", high = "#c53030", midpoint = 0,
      name = "NES"
    ) +
    scale_size_continuous(name = "Gene set size") +
    labs(title = title) +
    theme_embo(11) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    )
}

# -------------------- GO BP NES-centered PDF --------------------
p_bp_nes <- nes_centered_plot(gsea_go_bp, "GSEA GO BP (NES-centered: Up vs Down)")
save_pdf(p_bp_nes, "GSEA_GO_BP_NES_centered.pdf", w = 8.5, h = 7.5)

# -------------------- GO BP/CC/MF combined (NES-centered) --------------------
p_go_3panel_nes <- (
  nes_centered_plot(gsea_go_bp, "GO BP") |
    nes_centered_plot(gsea_go_cc, "GO CC") |
    nes_centered_plot(gsea_go_mf, "GO MF")
) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

save_pdf(p_go_3panel_nes, "GSEA_GO_BPCMF_NES_centered_3panel_EMBO.pdf", w = 16.5, h = 7.5)

# -------------------- KEGG + Reactome combined (NES-centered) --------------------
p_kegg <- nes_centered_plot(gsea_kegg,  "GSEA KEGG (NES-centered)")
p_react <- nes_centered_plot(gsea_react, "GSEA Reactome (NES-centered)")

p_kegg_react <- (p_kegg | p_react) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

save_pdf(p_kegg_react, "GSEA_KEGG_Reactome_NES_centered_2panel_EMBO.pdf", w = 15.5, h = 7.5)

# ==========================================================
# 2) Improved island/network plot
#    - subset top terms by p.adjust
#    - compute pairwise term similarity
#    - emapplot with minimal args (version safe)
#    - (optional) cluster labels via ggtitle + smaller text
# ==========================================================
make_emap_safe <- function(gsea_obj, title, keep_terms = network_keep_terms, layout = "kk") {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) return(NULL)
  
  tab <- as.data.frame(gsea_obj) |> arrange(p.adjust)
  keep_ids <- head(tab$ID, keep_terms)
  
  sub <- gsea_obj
  sub@result <- gsea_obj@result[gsea_obj@result$ID %in% keep_ids, , drop = FALSE]
  if (nrow(sub@result) < min_terms_for_network) return(NULL)
  
  sub_sim <- pairwise_termsim(sub)
  
  p <- tryCatch({
    emapplot(sub_sim, layout = layout, showCategory = nrow(sub_sim@result)) +
      ggtitle(title) +
      theme_embo(11) +
      theme(text = element_text(size = 8))
  }, error = function(e) {
    message("emapplot failed: ", conditionMessage(e))
    NULL
  })
  
  p
}

# GO BP network/islands (clearer by limiting to top terms)
p_bp_net <- make_emap_safe(gsea_go_bp, "GSEA GO BP – islands/network", keep_terms = 60, layout = "kk")
if (!is.null(p_bp_net)) {
  save_pdf(p_bp_net, "GSEA_GO_BP_islands_network.pdf", w = 10.5, h = 8.5)
} else {
  message("No GO BP network plot produced (not enough terms or emapplot failed).")
}

# KEGG / Reactome network plots (often fewer terms)
p_kegg_net <- make_emap_safe(gsea_kegg, "GSEA KEGG – islands/network", keep_terms = 50, layout = "kk")
if (!is.null(p_kegg_net)) save_pdf(p_kegg_net, "GSEA_KEGG_islands_network.pdf", w = 10.5, h = 8.5)

p_react_net <- make_emap_safe(gsea_react, "GSEA Reactome – islands/network", keep_terms = 50, layout = "kk")
if (!is.null(p_react_net)) save_pdf(p_react_net, "GSEA_Reactome_islands_network.pdf", w = 10.5, h = 8.5)

# -------------------- Done --------------------
message("\nDONE ✅ Enrichment outputs saved to:\n  ", outdir, "\n\nMain PDFs:\n",
        "  - GSEA_GO_BP_NES_centered.pdf\n",
        "  - GSEA_GO_BPCMF_NES_centered_3panel_EMBO.pdf\n",
        "  - GSEA_KEGG_Reactome_NES_centered_2panel_EMBO.pdf\n",
        "  - GSEA_GO_BP_islands_network.pdf (if enough terms)\n")
