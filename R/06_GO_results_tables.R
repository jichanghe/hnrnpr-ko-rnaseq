# ==========================================================
# GO + GSEA analysis from HNRNPR_EMBO.xlsx
# ==========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

# ---------------- INPUT ----------------
file <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/HNRNPR_EMBO.xlsx"
outdir <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/GO_results"
dir.create(outdir, showWarnings = FALSE)

# ---------------- READ FILE ----------------
df <- read_excel(file)
df <- as.data.frame(df)

# ---------------- CLEAN TABLE ----------------
df2 <- df %>%
  dplyr::select(gene, log2FoldChange, padj, pvalue) %>%
  filter(!is.na(gene), gene != "", !is.na(log2FoldChange)) %>%
  distinct(gene, .keep_all = TRUE)

# ---------------- MAP SYMBOL → ENTREZ ----------------
gene_map <- bitr(df2$gene,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)

df2 <- inner_join(df2, gene_map, by = c("gene" = "SYMBOL"))

cat("Mapped genes:", nrow(df2), "\n")

# ==========================================================
# 1️⃣ ORA (padj < 0.05)
# ==========================================================

sig_genes <- df2 %>%
  filter(padj < 0.05) %>%
  pull(ENTREZID)

ego_BP <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(as.data.frame(ego_BP),
          file.path(outdir, "ORA_GO_BP.csv"),
          row.names = FALSE)

pdf(file.path(outdir, "ORA_GO_BP_dotplot.pdf"), width=8, height=6)
dotplot(ego_BP, showCategory=20) +
  ggtitle("ORA GO Biological Process") +
  theme_bw()
dev.off()

# ==========================================================
# 2️⃣ GSEA (BEST METHOD)
# ==========================================================

geneList <- df2$log2FoldChange
names(geneList) <- df2$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

gsea_BP <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

write.csv(as.data.frame(gsea_BP),
          file.path(outdir, "GSEA_GO_BP.csv"),
          row.names = FALSE)

pdf(file.path(outdir, "GSEA_GO_BP_dotplot.pdf"), width=8, height=7)
dotplot(gsea_BP, showCategory=20) +
  ggtitle("GSEA GO Biological Process") +
  theme_bw()
dev.off()

cat("DONE ✓ GO results saved in:", outdir)
