# ===== PERFECT PCA plot (journal-grade) =====
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)

# --------- input / output ---------
file_path <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/pca/psi_pca.tsv"
out_pdf   <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/pca/psi_PCA_PC1_PC2_PERFECT.pdf"
out_png   <- "H:/JI3/HNRNPR/workdir_v8.1.0/results/pca/psi_PCA_PC1_PC2_PERFECT.png"

# --------- load PCA scores ---------
pca_df <- read_tsv(file_path, show_col_types = FALSE)

# first column = sample names
names(pca_df)[1] <- "Sample"

# numeric PCs
pca_df <- pca_df %>%
  mutate(across(starts_with("PC"), as.numeric))

# group assignment
pca_df <- pca_df %>%
  mutate(
    Group = case_when(
      grepl("^Ref", Sample, ignore.case = TRUE) ~ "Ref",
      grepl("^Alt", Sample, ignore.case = TRUE) ~ "Alt",
      TRUE ~ "Other"
    ),
    Group = factor(Group, levels = c("Ref", "Alt"))
  )

# --------- PERFECT PCA plot ---------
p <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  
  # main points
  geom_point(size = 5, alpha = 0.95) +
  
  # soft confidence ellipses (68% = 1 SD)
  stat_ellipse(
    aes(group = Group),
    type = "norm",
    level = 0.68,
    linewidth = 0.8,
    linetype = "solid",
    show.legend = FALSE
  ) +
  
  # sample labels
  geom_text(
    aes(label = Sample),
    vjust = -1.0,
    size = 4.2,
    fontface = "plain",
    show.legend = FALSE
  ) +
  
  # perfect color palette
  scale_color_manual(
    values = c(
      "Ref" = "#0072B2",
      "Alt" = "#D55E00"
    )
  ) +
  
  labs(
    title = "PCA of PSI Splicing Profiles",
    x = "PC1",
    y = "PC2"
  ) +
  
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title  = element_text(face = "bold"),
    axis.line   = element_line(linewidth = 0.9),
    axis.ticks  = element_line(linewidth = 0.9),
    legend.position = "none"
  )

# --------- save ---------
ggsave(out_pdf, p, width = 7, height = 6)
ggsave(out_png, p, width = 7, height = 6, dpi = 1000)

message("✅ PERFECT PCA saved to:\n  ", out_pdf, "\n  ", out_png)
