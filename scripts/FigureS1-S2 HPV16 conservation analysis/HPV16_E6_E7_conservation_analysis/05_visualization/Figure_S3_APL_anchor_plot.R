library(ggplot2)
library(readr)
library(dplyr)
library(stringr)

source("theme_publication.R")

in_dir  <- "../results/tables"
out_dir <- "../results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(file.path(in_dir, "APL_anchor_conservation.tsv"), show_col_types = FALSE)

df <- df %>%
  mutate(
    Mutation_label = paste0(APL, ": ", Mutation),
    Mutation_label = factor(
      Mutation_label,
      levels = c(
        "E6apl: T1Y",
        "E6apl: I2L",
        "E7apl1: L1Y",
        "E7apl2: T10V",
        "E7apl3: M2L",
        "E7apl3: T10V"
      )
    ),
    APL = factor(APL, levels = c("E6apl", "E7apl1", "E7apl2", "E7apl3"))
  )

apl_colors <- c(
  "E6apl"  = "#8B3A0F",
  "E7apl1" = "#E08A3C",
  "E7apl2" = "#C1622A",
  "E7apl3" = "#A64B1A"
)

p <- ggplot(df, aes(x = Mutation_label, y = WT_conservation_percent, fill = APL)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = paste0(sprintf("%.1f", WT_conservation_percent), "%")),
    vjust = -0.4, size = 4
  ) +
  scale_fill_manual(values = apl_colors, name = "APL") +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Conservation of WT residues targeted by APL substitutions",
    x = NULL,
    y = "WT residue conservation (%)"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

ggsave(
  filename = file.path(out_dir, "Figure_S3_APL_anchor_conservation.pdf"),
  plot = p, width = 8.2, height = 5.2
)

ggsave(
  filename = file.path(out_dir, "Figure_S3_APL_anchor_conservation.png"),
  plot = p, width = 8.2, height = 5.2, dpi = 600
)

print(p)