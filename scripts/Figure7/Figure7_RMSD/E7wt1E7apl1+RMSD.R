# ==========================================
# 结构免疫学 MD：50ns 完整版 RMSD 
# (极简纯净版：统计口径为 20-50ns，无背景阴影)
# ==========================================

library(ggplot2)
library(dplyr)

# 1. 设置工作目录 (请修改为你真实的绝对路径)
setwd("E:/HPV-T非临床文章/0.文章投稿+确定版本/Vaccines/三维结构图/HLA-A0201/图表制作/analysis/rmsd/data")

# 2. 读取文件 (无视原本错误的第一列)
wt_data <- read.table("E7wt1_9mer_rmsd.dat", col.names = c("Fake_Time", "RMSD"))
apl_data <- read.table("E7apl1_9mer_rmsd.dat", col.names = c("Fake_Time", "RMSD"))

# 3. 物理纠错：强制赋予真实的 50ns 物理时间 (0.005 ns/帧)
wt_data$Time_ns <- (1:nrow(wt_data)) * 0.005
apl_data$Time_ns <- (1:nrow(apl_data)) * 0.005

# 打上分组标签并合并数据集
wt_data$Group <- "WT (E7wt1)"
apl_data$Group <- "APL (E7apl1)"
merged_data <- bind_rows(wt_data, apl_data)

# --- 核心保留：平台期统计口径依然严格对齐 20-50 ns ---
plateau_stats <- merged_data %>%
  filter(Time_ns >= 20) %>%
  group_by(Group) %>%
  summarise(Mean = mean(RMSD), SD = sd(RMSD), .groups = "drop")

print(plateau_stats)
plateau_lab <- paste0(plateau_stats$Group, ":  ",
                      sprintf("%.2f", plateau_stats$Mean), " ± ",
                      sprintf("%.2f", plateau_stats$SD), " Å")

# 顶刊经典配色：WT 蓝 vs APL 橙
custom_colors <- c("WT (E7wt1)" = "#5BA3D0", "APL (E7apl1)" = "#E08A3C")

# 4. 开始高颜值绘图 (去除了阴影和虚线)
p <- ggplot(merged_data, aes(x = Time_ns, y = RMSD, color = Group)) +
  
  # 背景层：真实波动轨迹，调低透明度与线宽
  geom_line(alpha = 0.2, linewidth = 0.3) + 
  
  # 核心层：宏观收敛趋势拟合，加粗显示
  geom_smooth(method = "loess", span = 0.1, linewidth = 1.2, se = FALSE) +
  
  # 平台期均值 ± SD（带肽名+颜色，兼作图例，文案标明 20-50 ns）
  annotate("text", x = 1.5, y = 3.95, label = "Plateau RMSD (20-50 ns)",
           hjust = 0, vjust = 1, size = 3.6, fontface = "bold", color = "black") +
  annotate("text", x = 1.5, y = c(3.60, 3.28), label = plateau_lab,
           hjust = 0, vjust = 1, size = 4.0, fontface = "bold",
           color = custom_colors[plateau_stats$Group]) +
  
  # 强行画一根完美的 L 型黑线，取代系统默认的坐标轴
  annotate("path", x = c(0, 0, 50), y = c(4, 0, 0), color = "black", linewidth = 1.2, linejoin = "mitre") +
  
  theme_classic(base_size = 15) +
  scale_color_manual(values = custom_colors) +
  labs(x = "Time (ns)", y = "Backbone RMSD (Å)", color = "Peptide Type") +
  
  # 保留 0 刻度，并砍掉所有两端的延伸留白
  scale_x_continuous(breaks = seq(0, 50, by = 10), expand = c(0, 0)) + 
  scale_y_continuous(breaks = seq(0, 4, by = 1), expand = c(0, 0)) +
  
  # 允许在边界外绘图
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 4), clip = "off") + 
  
  # 顶刊排版细节调整
  theme(
    axis.line = element_blank(),
    axis.text.x = element_text(color = "black", size = 13, margin = margin(t = 6)),
    axis.text.y = element_text(color = "black", size = 13, margin = margin(r = 6)),
    axis.title = element_text(face = "bold", size = 16, margin = margin(t = 10, r = 10)),
    axis.ticks = element_line(linewidth = 1.0, color = "black"),
    axis.ticks.length = unit(0.2, "cm"), 
    legend.position = "none",
    legend.background = element_blank(),
    legend.key = element_blank(), 
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 12),
    plot.margin = margin(t = 15, r = 25, b = 15, l = 15, unit = "pt")
  )

# 在 RStudio 预览
print(p)

# 5. 导出为高清晰度 PDF 矢量图
out <- "Figure_RMSD_50ns_E7wt1_E7apl1_v1"
# PDF：矢量，cairo_pdf 能正确嵌入字体和 ± / Å 等字符
ggsave(paste0(out, ".pdf"), plot = p, width = 8.5, height = 6, device = cairo_pdf)

# PNG：600 dpi
ggsave(paste0(out, ".png"), plot = p, width = 8.5, height = 6, dpi = 600)

# TIFF：600 dpi + LZW 压缩（投稿主用）
ggsave(paste0(out, ".tif"), plot = p, width = 8.5, height = 6, dpi = 600, device = "tiff", compression = "lzw")