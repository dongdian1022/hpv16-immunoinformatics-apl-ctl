library(ggplot2)
library(patchwork)
library(dplyr)

# ══════════════════════════════════════════════════════════════════
# 数据
# ══════════════════════════════════════════════════════════════════
df <- data.frame(
  peptide = rep(c("E7wt1","E7wt2","E6wt",
                  "E7apl1","E7apl2","E7apl3","E6apl"), each = 3),
  round   = rep(c(1, 2, 3), times = 7),
  fold    = c(
    2.4,  9.1, 22.3,
    1.9, 16.3, 31.7,
    2.4, 15.4, 36.2,
    2.9,  9.3, 27.4,
    3.1, 14.6, 45.5,
    2.9, 15.7, 45.0,
    3.4, 15.7, 53.1
  ),
  group = rep(c("WT","WT","WT","APL","APL","APL","APL"), each = 3)
)

df$peptide <- factor(df$peptide,
                     levels = c("E7wt1","E7wt2","E6wt",
                                "E7apl1","E7apl2","E7apl3","E6apl"))
df$group <- factor(df$group, levels = c("WT","APL"))
df_r3    <- df %>% filter(round == 3)

# ══════════════════════════════════════════════════════════════════
# 配色
# ══════════════════════════════════════════════════════════════════
col_wt  <- "#3A7EBF"
col_apl <- "#C1622A"

pal <- c(
  "E7wt1"  = "#6AAED6",
  "E7wt2"  = "#3A7EBF",
  "E6wt"   = "#1A4F8A",
  "E7apl1" = "#E8A26A",
  "E7apl2" = "#C1622A",
  "E7apl3" = "#D4732A",
  "E6apl"  = "#8B3A0F"
)

# ══════════════════════════════════════════════════════════════════
# 共用 theme
# ══════════════════════════════════════════════════════════════════
vaccines_theme <- theme_classic(base_size = 12) +
  theme(
    axis.line          = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length  = unit(3.5, "pt"),
    axis.text.x        = element_text(size = 11, color = "#000000"),
    axis.text.y        = element_text(size = 11, color = "#000000"),
    axis.title.x       = element_text(size = 12, color = "#000000",
                                      face = "bold",
                                      margin = margin(t = 7)),
    axis.title.y       = element_text(size = 12, color = "#000000",
                                      face = "bold",
                                      margin = margin(r = 8)),
    panel.background   = element_rect(fill = "white", color = NA),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(8, 18, 8, 8),
    plot.tag           = element_text(size = 14, face = "bold",
                                      color = "#000000")
  )

# ══════════════════════════════════════════════════════════════════
# Fig 3A — 折线图
# ══════════════════════════════════════════════════════════════════
label_normal <- df %>%
  filter(round == 3,
         !peptide %in% c("E7apl2", "E7apl3"))

fold_apl2 <- 45.5
fold_apl3 <- 45.0

# 图例顺序：E6wt、E7wt1、E7wt2、E6apl、E7apl1、E7apl2、E7apl3
# 对应 group：WT、WT、WT、APL、APL、APL、APL
# linetype override 按此顺序：前3个dashed(WT)，后4个solid(APL)
legend_order   <- c("E6wt","E7wt1","E7wt2","E6apl","E7apl1","E7apl2","E7apl3")
legend_shape   <- c(21, 21, 21, 16, 16, 16, 16)      # WT空心，APL实心
legend_lty     <- c("dashed","dashed","dashed",
                    "solid","solid","solid","solid")  # WT虚线，APL实线

p_a <- ggplot(df, aes(x = round, y = fold,
                      color = peptide, group = peptide)) +
  
  geom_line(aes(linetype = group),
            linewidth = 1.0) +
  
  geom_point(aes(shape = group),
             size   = 2.4,
             stroke = 0.8,
             fill   = "white") +
  
  geom_text(data        = label_normal,
            aes(y       = fold, label = peptide),
            hjust       = -0.25,
            vjust       = 0.5,
            size        = 2.9,
            fontface    = "bold",
            show.legend = FALSE) +
  
  annotate("text",
           x = 3, y = fold_apl2 + 3.5,
           label = "E7apl2", hjust = -0.25, vjust = 0.5,
           size = 2.9, fontface = "bold", color = pal["E7apl2"]) +
  
  annotate("text",
           x = 3, y = fold_apl3 - 3.5,
           label = "E7apl3", hjust = -0.25, vjust = 0.5,
           size = 2.9, fontface = "bold", color = pal["E7apl3"]) +
  
  scale_color_manual(
    values = pal,
    name   = "Peptide",
    breaks = legend_order,
    guide  = guide_legend(
      order        = 1,
      override.aes = list(
        linewidth = 1.0,
        size      = 2.2,
        shape     = legend_shape,
        linetype  = legend_lty
      )
    )
  ) +
  
  scale_linetype_manual(values = c("WT" = "dashed", "APL" = "solid"),
                        guide  = "none") +
  scale_shape_manual(values   = c("WT" = 21, "APL" = 16),
                     guide    = "none") +
  
  scale_x_continuous(breaks = 1:3,
                     labels = c("Round 1", "Round 2", "Round 3"),
                     expand = expansion(mult = c(0.04, 0.30))) +
  
  scale_y_continuous(limits = c(0, 60),
                     breaks = seq(0, 60, by = 10),
                     expand = c(0, 0)) +
  
  labs(x = NULL, y = "Cumulative Fold Expansion", tag = "A") +
  
  vaccines_theme +
  theme(
    # 修1：图例移至左上，远离数据区；字号和间距紧凑
    legend.position   = c(0.16, 0.82),
    legend.title      = element_text(size = 9, face = "bold"),
    legend.text       = element_text(size = 8.5),
    legend.key.width  = unit(24, "pt"),   # 略缩小
    legend.key.height = unit(11, "pt"),   # 略缩小，整体更紧凑
    legend.background = element_rect(fill = "white", color = NA),
    legend.margin     = margin(3, 6, 3, 6)
  )

# ══════════════════════════════════════════════════════════════════
# Fig 3B — 分组柱状图
# ══════════════════════════════════════════════════════════════════
df_b <- data.frame(
  peptide = c("E7wt1","E7apl1",
              "E7wt2","E7apl2","E7apl3",
              "E6wt","E6apl"),
  fold    = c(22.3, 27.4,
              31.7, 45.5, 45.0,
              36.2, 53.1),
  group   = factor(c("WT","APL",
                     "WT","APL","APL",
                     "WT","APL"),
                   levels = c("WT","APL")),
  xpos    = c(1.0, 1.7,
              3.2, 3.9, 4.6,
              6.1, 6.8)
)

# 修3：分隔线颜色极淡
divider_x <- c(2.45, 5.35)

p_b <- ggplot(df_b, aes(x = xpos, y = fold, fill = group)) +
  
  geom_vline(xintercept = divider_x,
             linetype   = "dotted",
             color      = "#DDDDDD",     # 更淡
             linewidth  = 0.5) +
  
  geom_bar(stat  = "identity",
           width = 0.6,
           color = NA) +
  
  # 修5：柱顶数字间距加大
  geom_text(aes(y     = fold + 2.0,
                label = fold),
            size     = 3.2,
            fontface = "bold",
            color    = "#000000") +
  
  scale_fill_manual(values = c("WT" = col_wt, "APL" = col_apl),
                    name   = NULL) +
  
  scale_x_continuous(breaks = df_b$xpos,
                     labels = df_b$peptide,
                     expand = expansion(add = c(0.5, 0.5))) +
  
  scale_y_continuous(limits = c(0, 60),
                     breaks = seq(0, 60, by = 10),
                     expand = c(0, 0)) +
  
  labs(x = NULL, y = "Fold Expansion (Round 3)", tag = "B") +
  
  vaccines_theme +
  theme(
    axis.text.x       = element_text(size  = 10.5,
                                     color = "#000000",
                                     angle = 45,
                                     hjust = 1,
                                     vjust = 1),
    legend.position   = c(0.10, 0.91),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size   = unit(13, "pt"),
    legend.text       = element_text(size = 10),
    legend.margin     = margin(4, 8, 4, 8)
  )

# ══════════════════════════════════════════════════════════════════
# 拼图
# ══════════════════════════════════════════════════════════════════
fig3 <- p_a + p_b +
  plot_layout(widths = c(1.2, 1.0)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# ══════════════════════════════════════════════════════════════════
# 保存 — PNG / PDF / JPG / TIFF，全部 300 dpi
# ══════════════════════════════════════════════════════════════════
W <- 12   # 宽度 inches
H <- 5.5  # 高度 inches

# PNG — 300 dpi
ggsave("Figure3_AB_v9.png",
       plot   = fig3,
       width  = W, height = H,
       dpi    = 300,
       bg     = "white")

# PDF — 矢量（无 dpi 概念，打印质量最高）
ggsave("Figure3_AB_v9.pdf",
       plot   = fig3,
       width  = W, height = H,
       bg     = "white")

# JPG — 300 dpi，quality=95 减少压缩失真
ggsave("Figure3_AB_v9.jpg",
       plot    = fig3,
       width   = W, height = H,
       dpi     = 300,
       bg      = "white",
       device  = jpeg,
       quality = 95)

# TIFF — 300 dpi，LZW无损压缩（期刊投稿标准格式）
ggsave("Figure3_AB_v9.tiff",
       plot        = fig3,
       width       = W, height = H,
       dpi         = 300,
       bg          = "white",
       device      = tiff,
       compression = "lzw")

message("Done!  Figure3_AB_v9.png / .pdf / .jpg / .tiff")
print(fig3)

# install.packages(c("ggplot2", "patchwork", "dplyr"))