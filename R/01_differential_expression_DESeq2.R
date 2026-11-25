# ===== PERFECT Volcano Plot (DEG) =====
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

file_path <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/DEG.txt"
out_pdf   <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/DEG_volcano_x10.pdf"
out_png   <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/expression/DEG_volcano_x10.png"

padj_cut  <- 0.05
lfc_cut   <- 1
top_n_lab <- 10

df <- fread(file_path, data.table = FALSE, sep = "\t", header = TRUE) %>%
  mutate(
    gene_name = as.character(gene_name),
    log2FoldChange = as.numeric(log2FoldChange),
    padj = as.numeric(padj),
    neglog10_padj = -log10(pmax(padj, 1e-300)),
    group = case_when(
      padj < padj_cut & log2FoldChange >=  lfc_cut ~ "Up",
      padj < padj_cut & log2FoldChange <= -lfc_cut ~ "Down",
      TRUE                                         ~ "NS"
    ),
    group = factor(group, levels = c("Down","Up","NS"))
  )

label_df <- df %>%
  filter(group %in% c("Up","Down")) %>%
  group_by(group) %>%
  arrange(padj) %>%
  slice_head(n = top_n_lab) %>%
  ungroup()

p <- ggplot(df, aes(x = log2FoldChange, y = neglog10_padj)) +
  geom_point(aes(color = group), alpha = 0.85, size = 2.3) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_cut),     linetype = "dashed", linewidth = 0.4) +
  ggrepel::geom_label_repel(
    data = label_df,
    aes(label = gene_name, fill = group),
    color = "black",
    size = 3.5,
    label.size = 0.25,
    box.padding = 0.6,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Down" = "#0072B2",
    "Up"   = "#D55E00",
    "NS"   = "grey70"
  )) +
  scale_fill_manual(values = c(
    "Down" = "#A6CEE3",
    "Up"   = "#FDBF6F",
    "NS"   = "white"
  ), guide = "none") +
  coord_cartesian(xlim = c(-12, 10)) +   # 🔥 THIS IS THE KEY LINE
  labs(
    title = "Volcano Plot (DEGs)",
    x = "log2(Fold Change)",
    y = expression(-log[10](padj))
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

ggsave(out_pdf, p, width = 7.5, height = 6.5)
ggsave(out_png, p, width = 7.5, height = 6.5, dpi = 1000)
