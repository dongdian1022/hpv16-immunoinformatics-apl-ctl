# ==========================================
# 结构免疫学 MD 数据可视化：多肽 RMSF 柔性对比 (10-mer)
# (极简纯净版：完美对齐 20-50ns 数据与 9-mer 视觉风格)
# ==========================================

library(ggplot2)
library(dplyr)

# 1. 设置工作目录 (⚠️ 指向 rmsf_v3 文件夹)
setwd("E:/HPV-T非临床文章/0.文章投稿+确定版本/Vaccines/三维结构图/HLA-A0201/图表制作/analysis/rmsf/rmsf_v3")

# 2. 读取文件 (⚠️ 确保读取的是重算后的 20-50ns v3 纯净数据)
wt_data <- read.table("E7wt2_10mer_rmsf_20_50ns_v3.dat", col.names = c("Residue", "RMSF"))
apl_data <- read.table("E7apl2_10mer_rmsf_20_50ns_v3.dat", col.names = c("Residue", "RMSF"))

# 3. 数据预处理与位点对齐 (分配 P1-P10)
wt_data$Position <- 1:nrow(wt_data)
apl_data$Position <- 1:nrow(apl_data)

# 为数据打上分组标签
wt_data$Group <- "WT (E7wt2)"
apl_data$Group <- "APL (E7apl2)"

# 合并成主数据框
merged_data <- bind_rows(wt_data, apl_data)

# 4. 设定配色方案 (保持 10-mer 专属的蓝橙对比)
custom_colors <- c("WT (E7wt2)" = "#3A7EBF", "APL (E7apl2)" = "#C1622A")

# --- 锚位与改造位点标注参数 ---
anchor_pos <- c(2, 10)        # HLA-A*02:01 锚位 P2 与 C 端 P10
sub_pos    <- c(10)           # 改造位点
sub_lab    <- c("T10V")

# 5. 开始高颜值绘图
p <- ggplot(merged_data, aes(x = Position, y = RMSF, color = Group, group = Group)) +
  
  # 锚定口袋阴影 (置于底层)
  annotate("rect", xmin = anchor_pos - 0.45, xmax = anchor_pos + 0.45,
           ymin = 0, ymax = 3, fill = "grey80", alpha = 0.35) +
  
  # 改造位点：标签 + 向下三角指向该位点
  annotate("text", x = sub_pos, y = 2.80, label = sub_lab,
           fontface = "bold", size = 4.2, color = "black") +
  annotate("point", x = sub_pos, y = 2.65, shape = 25,
           fill = "black", color = "black", size = 2.6) +
  
  # 绘制点与线 
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5, shape = 16) + 
  
  # 🚨 魔法 1：复刻完美的 L 型黑线 (注意 X 轴坐标延伸到了 10.5)
  annotate("path", x = c(0.5, 0.5, 10.5), y = c(3, 0, 0), color = "black", linewidth = 1.2, linejoin = "mitre") +
  
  theme_classic(base_size = 15) +
  scale_color_manual(values = custom_colors) +
  
  # 坐标轴标签 
  labs(x = "Peptide Position", y = "Backbone RMSF (Å)") +
  
  # 🚨 魔法 2：严格控制 X 轴为 1-10，消除两端多余留白
  scale_x_continuous(breaks = 1:10, labels = paste0("P", 1:10), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 10.5), ylim = c(0, 3), clip = "off") + 
  
  # 顶刊排版细节调整
  theme(
    axis.line = element_blank(), # 废掉系统坐标轴，用我们画的 L 型线
    axis.text.x = element_text(color = "black", size = 13, face = "bold", margin = margin(t = 6)), 
    axis.text.y = element_text(color = "black", size = 13, margin = margin(r = 6)),
    axis.title = element_text(face = "bold", size = 16, margin = margin(t = 10, r = 10)),
    axis.ticks = element_line(linewidth = 1.0, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 图例悬浮于顶部居中，极致干净
    legend.position = c(0.5, 0.95),
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 12),
    plot.margin = margin(t = 15, r = 25, b = 15, l = 15, unit = "pt")
  )

# 预览图表
print(p)

# 6. 统一格式与尺寸导出
out <- "Figure_RMSF_20_50ns_E7wt2_E7apl2_v1"

ggsave(paste0(out, ".pdf"), plot = p, width = 8.5, height = 6)
ggsave(paste0(out, ".png"), plot = p, width = 8.5, height = 6, dpi = 600)
ggsave(paste0(out, ".tif"), plot = p, width = 8.5, height = 6, dpi = 600, device = "tiff", compression = "lzw")