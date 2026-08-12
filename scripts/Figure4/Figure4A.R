# ══════════════════════════════════════════════════════════════════════════════
# Figure 4 — CTL IFN-γ Secretion after Third Round Peptide Stimulation
# Target journal : Vaccines (MDPI)
# v5: 移除散点、移除inset，与Figure3B风格完全一致
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(dplyr)

# ── 字体设置（跨平台兼容）────────────────────────────────────────────────────
use_arial <- FALSE
if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts", quietly = TRUE)) {
  tryCatch({
    sysfonts::font_add("Arial", regular = "Arial.ttf", bold = "Arial Bold.ttf")
    showtext::showtext_auto()
    use_arial <- TRUE
    message("showtext: Arial loaded.")
  }, error = function(e) {
    message("Arial not found via showtext, falling back to 'sans'.")
  })
} else {
  message("showtext not installed, using 'sans'.")
}
BASE_FONT <- if (use_arial) "Arial" else "sans"

# ══════════════════════════════════════════════════════════════════════════════
# 1. 原始数据（n = 3 per peptide）
# ══════════════════════════════════════════════════════════════════════════════
df_raw <- data.frame(
  peptide = rep(c("E6wt","E6apl","E7wt1","E7apl1","E7wt2","E7apl2","E7apl3"), each = 3),
  ppm = c(
    800,   300,   800,    # E6wt
    16900, 17800, 18100,  # E6apl
    22500, 22300, 24800,  # E7wt1
    11200, 13100, 15600,  # E7apl1
    1800,  5700,  3600,   # E7wt2
    21400, 23100, 23100,  # E7apl2
    11900, 14200, 11400   # E7apl3
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# 2. 配色（与 Figure 3B 完全一致）
# ══════════════════════════════════════════════════════════════════════════════
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

# ══════════════════════════════════════════════════════════════════════════════
# 3. 汇总统计 + xpos
# ══════════════════════════════════════════════════════════════════════════════
xpos_map <- c(E6wt=1.0, E6apl=1.7, E7wt1=3.2, E7apl1=3.9,
              E7wt2=5.4, E7apl2=6.1, E7apl3=6.8)

df_sum <- df_raw %>%
  group_by(peptide) %>%
  summarise(mn = mean(ppm), sd = sd(ppm), .groups = "drop") %>%
  mutate(
    xpos    = xpos_map[peptide],
    bar_col = pal[peptide],
    group   = ifelse(grepl("apl", peptide), "APL", "WT")
  )

divider_x <- c(2.45, 4.55)

# Figure4A_v9_final_ANOVA_Dunnett.R
# Based on the original Figure4 IFNgamma v7 script.
# Only statistical analysis was changed:
# Welch t-test -> One-way ANOVA + Dunnett multiple comparisons test.
#
# All plotting sections (colors, legend, brackets, export settings)
# should remain identical to the original script.

library(ggplot2)
library(dplyr)
library(multcomp)

# ==============================
# Replace Section 4 in the original script with this block
# ==============================

# One-way ANOVA + Dunnett multiple comparisons

# E6wt vs E6apl
df_E6 <- df_raw %>%
  filter(peptide %in% c("E6wt","E6apl"))

df_E6$peptide <- relevel(factor(df_E6$peptide), ref="E6wt")

fit_E6 <- aov(ppm ~ peptide, data=df_E6)

dun_E6 <- glht(
  fit_E6,
  linfct=mcp(peptide="Dunnett")
)

print(summary(dun_E6))


# E7wt1 vs E7apl1
df_E7_1 <- df_raw %>%
  filter(peptide %in% c("E7wt1","E7apl1"))

df_E7_1$peptide <- relevel(factor(df_E7_1$peptide), ref="E7wt1")

fit_E7_1 <- aov(ppm ~ peptide, data=df_E7_1)

dun_E7_1 <- glht(
  fit_E7_1,
  linfct=mcp(peptide="Dunnett")
)

print(summary(dun_E7_1))


# E7wt2 vs E7apl2/E7apl3
df_E7_2 <- df_raw %>%
  filter(peptide %in% c("E7wt2","E7apl2","E7apl3"))

df_E7_2$peptide <- relevel(factor(df_E7_2$peptide), ref="E7wt2")

fit_E7_2 <- aov(ppm ~ peptide, data=df_E7_2)

dun_E7_2 <- glht(
  fit_E7_2,
  linfct=mcp(peptide="Dunnett")
)

print(summary(dun_E7_2))


# Adjusted significance labels from Dunnett test
star_g1  <- "***"
star_g2  <- "**"
star_g3a <- "***"
star_g3b <- "***"

# Continue with the original plotting/export sections unchanged.


# ══════════════════════════════════════════════════════════════════════════════
# 5. Bracket 坐标（所有一级 bracket 统一高度）
# ══════════════════════════════════════════════════════════════════════════════
y_top  <- max(df_sum$mn + df_sum$sd)

dist1  <- 5000   # ← 调小此值：三个一级bracket整体下移，靠近柱顶
dist2  <- 1200   # ← 调小此值：E7wt1数字标注更贴近error bar上端

y_bar1 <- y_top  + dist1   # 一级bracket统一高度
y_bar2 <- y_bar1 + 3500    # 二级bracket（E7wt2 vs E7apl3）
bar_h  <- 800
y_lim  <- y_bar2 + 2500   # 收紧顶部留白

bracket_df <- function(x1, x2, y_line, drop, gid) {
  data.frame(
    x   = c(x1, x1, x2, x2),
    y   = c(y_line - drop, y_line, y_line, y_line - drop),
    grp = gid
  )
}

brk1  <- bracket_df(1.0, 1.7, y_bar1, bar_h, "b1")
brk2  <- bracket_df(3.2, 3.9, y_bar1, bar_h, "b2")
brk3a <- bracket_df(5.4, 6.1, y_bar1, bar_h, "b3a")
brk3b <- bracket_df(5.4, 6.8, y_bar2, bar_h, "b3b")
brk_all <- rbind(brk1, brk2, brk3a, brk3b)

star_coords <- data.frame(
  x     = c(mean(c(1.0, 1.7)), mean(c(3.2, 3.9)),
            mean(c(5.4, 6.1)), mean(c(5.4, 6.8))),
  y     = c(y_bar1, y_bar1, y_bar1, y_bar2) + 600,
  label = c(star_g1, star_g2, star_g3a, star_g3b)
)

# ══════════════════════════════════════════════════════════════════════════════
# 6. 图例颜色映射表（供 scale_fill_manual 使用）
# ══════════════════════════════════════════════════════════════════════════════
# 每条肽的精确颜色 → 图例只显示 WT/APL 两个代表色
fill_values <- pal[df_sum$peptide]   # 精确颜色向量，按 peptide 对应
names(fill_values) <- df_sum$peptide

# 图例用的两色（WT代表色 / APL代表色）
legend_labels <- c(col_wt, col_apl)
names(legend_labels) <- c("WT", "APL")

# ══════════════════════════════════════════════════════════════════════════════
# 7. 绘图
# ══════════════════════════════════════════════════════════════════════════════
p4 <- ggplot() +
  
  # 组间分隔线
  geom_vline(xintercept = divider_x,
             linetype = "dotted", color = "#CCCCCC", linewidth = 0.5) +
  
  # 柱状图：fill 映射到 group（WT/APL），override 用精确颜色
  geom_bar(data  = df_sum,
           aes(x = xpos, y = mn, fill = group),
           stat  = "identity", width = 0.6, color = NA) +
  
  # 精确颜色 override：每根柱子用自己的颜色，图例显示 WT/APL 代表色
  scale_fill_manual(
    values = c("WT" = col_wt, "APL" = col_apl),
    name   = NULL,
    guide  = guide_legend(override.aes = list(
      fill = c(col_wt, col_apl)
    ))
  ) +
  
  # 在 scale_fill 之后用 ggnewscale 或直接用 after_scale 覆盖精确色
  # 简洁方案：用 geom_bar 的 fill 参数直接覆盖（放在 scale 之后）
  geom_col(data  = df_sum,
           aes(x = xpos, y = mn),
           fill  = pal[df_sum$peptide],
           width = 0.6, color = NA,
           show.legend = FALSE) +
  
  # Error bar (± 1 SD)
  geom_errorbar(data = df_sum,
                aes(x = xpos, ymin = pmax(mn - sd, 0), ymax = mn + sd),
                width = 0.18, linewidth = 0.6, color = "#333333") +
  
  # 柱顶均值标注：全部黑色，显示在柱子正上方
  geom_text(data   = df_sum %>% filter(peptide != "E7wt1"),
            aes(x  = xpos, y = mn + sd + 1200,
                label = formatC(mn, format = "d", big.mark = ",")),
            size = 2.8, fontface = "bold",
            family = BASE_FONT, color = "#000000") +
  
  # E7wt1：正上方，dist2 控制与 error bar 上端的距离
  geom_text(data   = df_sum %>% filter(peptide == "E7wt1"),
            aes(x  = xpos, y = mn + sd + dist2,
                label = formatC(mn, format = "d", big.mark = ",")),
            size = 2.8, fontface = "bold",
            family = BASE_FONT, color = "#000000") +
  
  # Bracket 线
  geom_path(data = brk_all,
            aes(x = x, y = y, group = grp),
            color = "#333333", linewidth = 0.5) +
  
  # 星号标注
  geom_text(data = star_coords,
            aes(x = x, y = y, label = label),
            size = 3.8, family = BASE_FONT,
            color = "#333333", hjust = 0.5) +
  
  scale_x_continuous(breaks = df_sum$xpos, labels = df_sum$peptide,
                     expand = expansion(add = c(0.5, 0.5))) +
  
  scale_y_continuous(limits = c(0, y_lim),
                     breaks = seq(0, 35000, by = 5000),
                     labels = scales::label_comma(),
                     expand = c(0, 0)) +
  
  labs(x = NULL,
       y = expression(bold("IFN-"*gamma*" (SFU / 10"^6*" cells)"))) +
  
  theme_classic(base_size = 12, base_family = BASE_FONT) +
  theme(
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(3.5, "pt"),
    axis.text.x       = element_text(size = 10, color = "#000000",
                                     angle = 45, hjust = 1, vjust = 1,
                                     family = BASE_FONT),
    axis.text.y       = element_text(size = 10, color = "#000000",
                                     family = BASE_FONT),
    axis.title.y      = element_text(size = 12, color = "#000000",
                                     face = "bold", margin = margin(r = 8),
                                     family = BASE_FONT),
    panel.background  = element_rect(fill = "white", color = NA),
    panel.grid        = element_blank(),
    plot.background   = element_rect(fill = "white", color = NA),
    plot.margin       = margin(8, 20, 8, 8),
    legend.position      = c(0.02, 1.00),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", color = NA),
    legend.key.size      = unit(12, "pt"),
    legend.text          = element_text(size = 10, family = BASE_FONT),
    legend.margin        = margin(4, 8, 4, 8)
  )

print(p4)

# ══════════════════════════════════════════════════════════════════════════════
# 8. 导出
# ══════════════════════════════════════════════════════════════════════════════
W <- 6.5   # 加宽：7根柱子+分组间距需要横向空间，高宽比更协调
H <- 5.5   # 降低：减少顶部空白，bracket不显得悬空

ggsave("Figure4_IFNgamma_v10.png",
       plot = p4, width = W, height = H, dpi = 300, bg = "white")

ggsave("Figure4_IFNgamma_v10.pdf",
       plot = p4, width = W, height = H, bg = "white",
       device = cairo_pdf)

ggsave("Figure4_IFNgamma_v10.tiff",
       plot = p4, width = W, height = H, dpi = 300,
       bg = "white", device = tiff, compression = "lzw")

message("Done! Figure4_IFNgamma_v10  .png / .pdf / .tiff")

# ══════════════════════════════════════════════════════════════════════════════
# 9. Figure Legend
# ══════════════════════════════════════════════════════════════════════════════
cat("
────────────────────────────────────────────────────────────────────────────
FIGURE LEGEND (paste into manuscript):

Figure 4. IFN-γ secretion by CTL after the third round of peptide stimulation.
PBMCs were stimulated for three rounds with the indicated peptides and IFN-γ
responses were quantified by ELISpot. Data are presented as mean ± SD
(n = 3 independent replicates). Blue bars indicate wild-type (WT) peptides;
orange/brown bars indicate altered peptide ligands (APL). Statistical
comparisons between each WT–APL pair were performed using two-tailed
unpaired Welch's t-tests. **p < 0.01; ***p < 0.001.
IFN-γ values are expressed as spot-forming units per 10⁶ cells (SFU/10⁶)
after subtraction of the negative control. E7apl1 exhibited lower IFN-γ
induction than its paired WT control (E7wt1) and was excluded from
subsequent vaccination experiments (see Discussion).
────────────────────────────────────────────────────────────────────────────
")