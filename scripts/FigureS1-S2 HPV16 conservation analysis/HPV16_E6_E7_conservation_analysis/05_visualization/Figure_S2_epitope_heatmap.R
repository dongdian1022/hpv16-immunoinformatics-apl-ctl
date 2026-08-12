library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)

source("theme_publication.R")

in_dir  <- "../results/tables"
out_dir <- "../results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(file.path(in_dir, "epitope_conservation.tsv"), show_col_types = FALSE)

# 统一排序
peptide_order <- c("E6wt", "E6apl", "E7wt1", "E7apl1", "E7wt2", "E7apl2", "E7apl3")

df <- df %>%
  mutate(
    Peptide = factor(Peptide, levels = peptide_order),
    Residue_position = as.integer(Residue_position)
  )

# 补齐1~10位，9-mer的P10补NA，便于画整齐热图
all_grid <- expand.grid(
  Peptide = peptide_order,
  Residue_position = 1:10,
  stringsAsFactors = FALSE
)

df_plot <- all_grid %>%
  left_join(df, by = c("Peptide", "Residue_position")) %>%
  mutate(
    Peptide = factor(Peptide, levels = rev(peptide_order)),
    Pos_label = paste0("P", Residue_position),
    cell_label = ifelse(
      is.na(Conservation_percent),
      "",
      paste0(Reference_AA, "\n", sprintf("%.1f", Conservation_percent))
    )
  )

p <- ggplot(df_plot, aes(x = factor(Residue_position), y = Peptide, fill = Conservation_percent)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = cell_label), size = 3.1, lineheight = 0.9) +
  scale_fill_gradient(
    low = "#F5D0C5",
    high = "#1A4F8A",
    na.value = "grey92",
    limits = c(70, 100),
    name = "Conservation (%)"
  ) +
  scale_x_discrete(labels = paste0("P", 1:10)) +
  labs(
    title = "Epitope-level conservation of HPV16 E6/E7 WT epitopes and APL candidates",
    x = "Peptide residue position",
    y = NULL
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "Figure_S2_epitope_heatmap.pdf"),
  plot = p, width = 9, height = 5.8
)

ggsave(
  filename = file.path(out_dir, "Figure_S2_epitope_heatmap.png"),
  plot = p, width = 9, height = 5.8, dpi = 600
)

print(p)