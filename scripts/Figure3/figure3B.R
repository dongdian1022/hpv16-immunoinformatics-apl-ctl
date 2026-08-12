library(ggplot2)

# ══════════════════════════════════════════════════════════════════
# 原始数据（n=3，保留一位小数）
# ══════════════════════════════════════════════════════════════════
df_raw <- data.frame(
  peptide = rep(c("E7wt2","E7apl2","E7apl3",
                  "E6wt","E6apl",
                  "E7wt1","E7apl1"), each = 3),
  fold = round(c(
    27.92, 32.65, 34.53,
    41.25, 45.88, 49.37,
    40.62, 45.91, 48.47,
    32.06, 36.94, 39.60,
    47.23, 54.69, 57.38,
    19.85, 22.74, 24.31,
    24.16, 27.83, 30.21
  ), 1)
)

# ══════════════════════════════════════════════════════════════════
# 配色
# ══════════════════════════════════════════════════════════════════
pal <- c(
  "E7wt1"  = "#6AAED6",
  "E7wt2"  = "#3A7EBF",
  "E6wt"   = "#1A4F8A",
  "E7apl1" = "#E8A26A",
  "E7apl2" = "#C1622A",
  "E7apl3" = "#D4732A",
  "E6apl"  = "#8B3A0F"
)

col_wt  <- "#3A7EBF"
col_apl <- "#C1622A"

# ══════════════════════════════════════════════════════════════════
# 手动建完整数据框
# 顺序：E6wt-E6apl / E7wt1-E7apl1 / E7wt2-E7apl2-E7apl3
# ══════════════════════════════════════════════════════════════════
calc_mean <- function(pep) round(mean(df_raw$fold[df_raw$peptide == pep]), 1)
calc_sd   <- function(pep) round(sd(df_raw$fold[df_raw$peptide == pep]),   1)

pep_order <- c("E6wt","E6apl","E7wt1","E7apl1","E7wt2","E7apl2","E7apl3")

df_b <- data.frame(
  peptide = pep_order,
  group   = factor(
    c("WT","APL","WT","APL","WT","APL","APL"),
    levels = c("WT","APL")
  ),
  xpos    = c(1.0, 1.7, 3.2, 3.9, 5.4, 6.1, 6.8),
  mean    = sapply(pep_order, calc_mean),
  sd      = sapply(pep_order, calc_sd),
  bar_col = unname(pal[pep_order])   # 修2：确保颜色与肽名精确对应
)

print(df_b[, c("peptide","mean","sd","bar_col")])

divider_x <- c(2.45, 4.55)

# ══════════════════════════════════════════════════════════════════
# 统计检验
# ══════════════════════════════════════════════════════════════════
get_star <- function(pep1, pep2) {
  x1 <- df_raw$fold[df_raw$peptide == pep1]
  x2 <- df_raw$fold[df_raw$peptide == pep2]
  p  <- t.test(x1, x2, paired = FALSE)$p.value
  if      (p < 0.001) "***"
  else if (p < 0.01)  "**"
  else if (p < 0.05)  "*"
  else                "ns"
}

star_g1  <- get_star("E6wt",  "E6apl")
star_g3a <- get_star("E7wt2", "E7apl2")
star_g3b <- get_star("E7wt2", "E7apl3")

cat("E6wt  vs E6apl: ", star_g1,  "\n")
cat("E7wt2 vs E7apl2:", star_g3a, "\n")
cat("E7wt2 vs E7apl3:", star_g3b, "\n")

# ══════════════════════════════════════════════════════════════════
# 统计标注 Y 坐标
# 修4：Y轴上限动态收紧
# ══════════════════════════════════════════════════════════════════
y_top  <- max(df_b$mean + df_b$sd)
y_bar1 <- y_top  + 3.0
y_bar2 <- y_bar1 + 4.0
bar_h  <- 0.7
y_lim  <- y_bar2 + 3.5   # 修4：收紧上限，减少顶部留白

# ══════════════════════════════════════════════════════════════════
# bracket 辅助函数
# ══════════════════════════════════════════════════════════════════
bracket_df <- function(x1, x2, y_line, drop, group_id) {
  data.frame(
    x   = c(x1, x1, x2, x2),
    y   = c(y_line - drop, y_line, y_line, y_line - drop),
    grp = group_id
  )
}

brk1  <- bracket_df(1.0, 1.7, y_bar1, bar_h, "b1")
brk3a <- bracket_df(5.4, 6.1, y_bar1, bar_h, "b3a")
brk3b <- bracket_df(5.4, 6.8, y_bar2, bar_h, "b3b")
brk_all <- rbind(brk1, brk3a, brk3b)

# ══════════════════════════════════════════════════════════════════
# 图例数据（手动构造，保证WT/APL图例正常显示）
# 修1：用独立的 legend_df 驱动图例，与柱子颜色解耦
# ══════════════════════════════════════════════════════════════════
legend_df <- data.frame(
  group = factor(c("WT","APL"), levels = c("WT","APL")),
  x     = c(-99, -99),   # 画到图外，不可见
  y     = c(-99, -99)
)

# ══════════════════════════════════════════════════════════════════
# 画图
# ══════════════════════════════════════════════════════════════════
p_b <- ggplot() +
  
  # 组间分隔线
  geom_vline(xintercept = divider_x,
             linetype   = "dotted",
             color      = "#CCCCCC",
             linewidth  = 0.5) +
  
  # 柱状图：每根柱子单独指定颜色
  geom_bar(data  = df_b,
           aes(x = xpos, y = mean),
           stat  = "identity",
           width = 0.6,
           fill  = df_b$bar_col,
           color = NA) +
  
  # Error bar
  geom_errorbar(data = df_b,
                aes(x    = xpos,
                    ymin = mean - sd,
                    ymax = mean + sd),
                width     = 0.18,
                linewidth = 0.6,
                color     = "#333333") +
  
  # 柱顶均值，修3：字号从3.3缩小到3.0
  geom_text(data = df_b,
            aes(x     = xpos,
                y     = mean + sd + 1.4,
                label = sprintf("%.1f", mean)),
            size     = 3.0,
            fontface = "bold",
            color    = "#000000") +
  
  # 无缝 bracket
  geom_path(data        = brk_all,
            aes(x = x, y = y, group = grp),
            color       = "#333333",
            linewidth   = 0.5) +
  
  # 星号标注
  annotate("text", x = 1.35, y = y_bar1 + 0.7,
           label = star_g1,  size = 3.8, color = "#333333", hjust = 0.5) +
  annotate("text", x = 5.75, y = y_bar1 + 0.7,
           label = star_g3a, size = 3.8, color = "#333333", hjust = 0.5) +
  annotate("text", x = 6.10, y = y_bar2 + 0.7,
           label = star_g3b, size = 3.8, color = "#333333", hjust = 0.5) +
  
  # 修1：用不可见柱驱动 WT/APL 图例
  geom_bar(data  = legend_df,
           aes(x = x, y = y, fill = group),
           stat  = "identity") +
  scale_fill_manual(values = c("WT" = col_wt, "APL" = col_apl),
                    name   = NULL) +
  
  scale_x_continuous(breaks = df_b$xpos,
                     labels = df_b$peptide,
                     expand = expansion(add = c(0.5, 0.5))) +
  
  scale_y_continuous(limits = c(0, y_lim),
                     breaks = seq(0, 60, by = 10),
                     expand = c(0, 0)) +
  
  labs(x = NULL,
       y = "Fold Expansion (Round 3)") +
  
  theme_classic(base_size = 12) +
  theme(
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(3.5, "pt"),
    axis.text.x       = element_text(size  = 10.5, color = "#000000",
                                     angle = 45, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size  = 11, color = "#000000"),
    axis.title.y      = element_text(size  = 12, color = "#000000",
                                     face  = "bold", margin = margin(r = 8)),
    panel.background  = element_rect(fill = "white", color = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.background   = element_rect(fill = "white", color = NA),
    plot.margin       = margin(8, 18, 8, 8),
    legend.position      = c(0.02, 1.00),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", color = NA),
    legend.key.size      = unit(13, "pt"),
    legend.text          = element_text(size = 10),
    legend.margin        = margin(4, 8, 4, 8)
  )

print(p_b)

# ══════════════════════════════════════════════════════════════════
# 保存
# ══════════════════════════════════════════════════════════════════
W <- 5.5
H <- 6.2

ggsave("Figure3B_final-v3.png",  plot = p_b, width = W, height = H,
       dpi = 300, bg = "white")

ggsave("Figure3B_final-v3.pdf",  plot = p_b, width = W, height = H,
       bg = "white")

ggsave("Figure3B_final-v3.jpg",  plot = p_b, width = W, height = H,
       dpi = 300, bg = "white", device = jpeg, quality = 95)

ggsave("Figure3B_final-v3.tiff", plot = p_b, width = W, height = H,
       dpi = 300, bg = "white", device = tiff, compression = "lzw")

message("Done!  Figure3B_final.png / .pdf / .jpg / .tiff")