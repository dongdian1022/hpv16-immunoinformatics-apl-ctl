library(ggplot2)
library(readr)
library(dplyr)
library(scales)

source("theme_publication.R")

in_dir  <- "../results/tables"
out_dir <- "../results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

e6 <- read_tsv(file.path(in_dir, "E6_full_length_conservation.tsv"), show_col_types = FALSE)
e7 <- read_tsv(file.path(in_dir, "E7_full_length_conservation.tsv"), show_col_types = FALSE)

# 统一列名兼容
colnames(e6) <- c("Gene", "Position", "Major_AA", "Conservation_percent", "N_sequences")
colnames(e7) <- c("Gene", "Position", "Major_AA", "Conservation_percent", "N_sequences")

e6$Gene <- "E6"
e7$Gene <- "E7"

df <- bind_rows(e6, e7)

highlight_df <- tibble::tribble(
  ~Gene, ~xmin, ~xmax, ~label, ~fill_col, ~line_col,
  "E6",  29,   38,   "E6wt / E6apl\n(29–38)", "#1A4F8A", "#1A4F8A",
  "E7",  11,   20,   "E7wt2 / E7apl2 / E7apl3\n(11–20)", "#3A7EBF", "#3A7EBF",
  "E7",  82,   90,   "E7wt1 / E7apl1\n(82–90)", "#5BA3D0", "#5BA3D0"
)

p <- ggplot(df, aes(x = Position, y = Conservation_percent)) +
  geom_rect(
    data = highlight_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = alpha(highlight_df$fill_col, 0.15),
    color = NA
  ) +
  geom_line(linewidth = 0.7, colour = "black") +
  geom_hline(yintercept = seq(70, 100, 5), colour = "grey90", linewidth = 0.3) +
  facet_wrap(~Gene, ncol = 1, scales = "free_x") +
  scale_y_continuous(
    limits = c(70, 100),
    breaks = seq(70, 100, 5),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Full-length conservation of HPV16 E6 and E7",
    x = "Amino acid position",
    y = "Conservation (%)"
  ) +
  theme_publication() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

# 添加文字标注
p <- p +
  annotate("text", x = 33.5, y = 71.5, label = "E6wt / E6apl\n(29–38)",
           colour = "#1A4F8A", size = 4, fontface = "bold") +
  annotate("text", x = 15.5, y = 71.5, label = "E7wt2 / E7apl2 / E7apl3\n(11–20)",
           colour = "#3A7EBF", size = 4, fontface = "bold") +
  annotate("text", x = 86, y = 71.5, label = "E7wt1 / E7apl1\n(82–90)",
           colour = "#5BA3D0", size = 4, fontface = "bold")

ggsave(
  filename = file.path(out_dir, "Figure_S1_full_length_conservation.pdf"),
  plot = p, width = 8.5, height = 8.5
)

ggsave(
  filename = file.path(out_dir, "Figure_S1_full_length_conservation.png"),
  plot = p, width = 8.5, height = 8.5, dpi = 600
)

print(p)