# ==========================================
# 结构免疫学 MD：多肽-HLA 氢键全自动精准捕获与绘图
# 【Nature/JITC 顶刊终极发表版】
# 特性：内部氢键过滤 + 0值补齐排版 + 高级标签清洗 + 动态虚线分区
# 针对 10-mer 体系绝对编号：P1(376), P2(377), P9(385)
# 针对20-50 ns的窗口
# ==========================================

library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

# ---------------------------------------------------------
# 1. 读取原始数据 (防弹读取模式)
# ---------------------------------------------------------
# 🚨 请确保修改为你电脑上的真实路径
setwd("E:/HPV-T非临床文章/0.文章投稿+确定版本/Vaccines/三维结构图/HLA-A0201/图表制作/analysis/hbond/data-v2")

wt_raw <- read.table("E7wt1_9mer_hbond_avg_20_50ns.dat", header = TRUE, comment.char = "", check.names = FALSE)
apl_raw <- read.table("E7apl1_9mer_hbond_avg_20_50ns.dat", header = TRUE, comment.char = "", check.names = FALSE)

# 清理带有 # 的表头
colnames(wt_raw)[1] <- sub("^#\\s*", "", colnames(wt_raw)[1])
colnames(apl_raw)[1] <- sub("^#\\s*", "", colnames(apl_raw)[1])

# 赋予分组标签
wt_raw$Group <- "WT (E7wt1)"
apl_raw$Group <- "APL (E7apl1)"
all_raw <- bind_rows(wt_raw, apl_raw)

# ---------------------------------------------------------
# 2. 核心捞取逻辑与高级标签清洗
# ---------------------------------------------------------
target_p1 <- "376"
target_p2 <- "377"
target_p9 <- "384"

clean_data <- all_raw %>%
  filter(
    grepl(paste0("_", target_p1, "@"), Acceptor) | grepl(paste0("_", target_p1, "@"), Donor) |
      grepl(paste0("_", target_p2, "@"), Acceptor) | grepl(paste0("_", target_p2, "@"), Donor) |
      grepl(paste0("_", target_p9, "@"), Acceptor) | grepl(paste0("_", target_p9, "@"), Donor)
  ) %>%
  rowwise() %>%
  mutate(
    Acc_Res = sub("@.*", "", Acceptor),
    Don_Res = sub("@.*", "", Donor),
    
    Peptide_Pos = case_when(
      grepl(paste0("_", target_p1, "$"), Acc_Res) | grepl(paste0("_", target_p1, "$"), Don_Res) ~ "P1",
      grepl(paste0("_", target_p2, "$"), Acc_Res) | grepl(paste0("_", target_p2, "$"), Don_Res) ~ "P2",
      grepl(paste0("_", target_p9, "$"), Acc_Res) | grepl(paste0("_", target_p9, "$"), Don_Res) ~ "P9",
      TRUE ~ "Unknown"
    ),
    
    Other_Res = ifelse(grepl(paste0("_(", target_p1, "|", target_p2, "|", target_p9, ")$"), Acc_Res), Don_Res, Acc_Res),
    Other_Num = as.numeric(str_extract(Other_Res, "\\d+")),
    
    # 🚨 标签整容术：把生硬的 TYR_159 变成高大上的 Tyr159
    Clean_AA = paste0(str_to_title(str_extract(Other_Res, "[A-Za-z]+")), Other_Num),
    
    # 生成分为上下两行、极具呼吸感的高级标签
    Plot_Label = paste0(Clean_AA, "\n(", Peptide_Pos, ")")
  ) %>%
  ungroup() %>%
  # 彻底过滤多肽内部氢键和水分子
  filter(Other_Num < 300)

# ---------------------------------------------------------
# 3. 数据合并与 0 值强制补齐排版
# ---------------------------------------------------------
final_plot_data <- clean_data %>%
  group_by(Group, Peptide_Pos, Plot_Label) %>%
  summarise(Occupancy = min(sum(Frac * 100), 100), .groups = "drop") %>%
  filter(Occupancy >= 10) %>%
  mutate(Group = factor(Group, levels = c("WT (E7wt1)", "APL (E7apl1)"))) %>%
  # 强制补齐 0 值，确保红灰柱子同等纤细
  complete(Group, nesting(Peptide_Pos, Plot_Label), fill = list(Occupancy = 0)) %>%
  mutate(Peptide_Pos = factor(Peptide_Pos, levels = c("P1", "P2", "P9"))) %>%
  arrange(Peptide_Pos, desc(Occupancy))

# 锁定 X 轴绘图顺序
final_plot_data$Plot_Label <- factor(final_plot_data$Plot_Label, levels = unique(final_plot_data$Plot_Label))

write.csv(final_plot_data, "Hbond_Final_E7WT1_E7APL1.csv", row.names = FALSE)
cat("\n✅ 终极数据处理完毕！请准备迎接绝美图表出炉！\n\n")

# ---------------------------------------------------------
# 4. 终极顶刊绘图 (极简高级感 + 科学严谨性修正)
# ---------------------------------------------------------
# 动态计算 P1 和 P2 分界线的位置
n_P1 <- sum(grepl("\\(P1\\)", levels(final_plot_data$Plot_Label)))
n_P2 <- sum(grepl("\\(P2\\)", levels(final_plot_data$Plot_Label)))

custom_colors <- c("WT (E7wt1)" = "#5BA3D0", "APL (E7apl1)" = "#E08A3C")

p <- ggplot(final_plot_data, aes(x = Plot_Label, y = Occupancy, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7, color = "black", linewidth = 0.3) +
  
  geom_vline(xintercept = n_P1 + 0.5, linetype = "dashed", color = "grey75", linewidth = 0.6) +
  geom_vline(xintercept = n_P1 + n_P2 + 0.5, linetype = "dashed", color = "grey75", linewidth = 0.6) +
  
  # 🚨 修正 1：L型直角黑线的 Y 轴最高点死死锁定在 100！
  annotate("path", x = c(0.4, 0.4, length(levels(final_plot_data$Plot_Label)) + 0.6), y = c(100, 0, 0), 
           color = "black", linewidth = 1.2, linejoin = "mitre") +
  
  theme_classic(base_size = 15) +
  scale_fill_manual(values = custom_colors) +
  
  # 🚨 修正 2：删除累赘的 X 轴标题，保持极简高级感
  labs(x = NULL, y = "H-bond Occupancy (%)") +
  
  # 🚨 修正 3：Y 轴视口严格限制在 0 到 100 的物理极限
  scale_y_continuous(breaks = seq(0, 100, by = 20), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100), clip = "off") + 
  
  theme(
    axis.line = element_blank(),
    
    # 因为删除了 X 轴标题，这里稍微增加 margin 让标签和底边保持优雅距离
    axis.text.x = element_text(color = "black", size = 12, face = "bold", hjust = 0.5, margin = margin(t = 10)),
    
    axis.text.y = element_text(color = "black", size = 13, margin = margin(r = 6)),
    
    # 删除了 X 轴标题，所以也不需要设置 axis.title.x 了
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 15)),
    
    axis.ticks = element_line(linewidth = 1.0, color = "black"),
    axis.ticks.x = element_blank(), 
    axis.ticks.length.y = unit(0.2, "cm"), 
    
    legend.position = c(0.25, 1.00),   # 往左减小x、往上增大y(最高1.0)；如挡住柱子就回调
    legend.justification = c(0.5, 1),
    legend.direction = "vertical",
    legend.margin = margin(0, 0, 0, 0),     # 去掉图例自带内边距，让它紧贴顶部
    legend.key.spacing.y = unit(10, "pt"),   # 上下两个色块之间的间距，想更远就调大
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 10),
    plot.caption = element_text(hjust = 0, size = 9),
    legend.background = element_blank(),
    
    plot.margin = margin(t = 20, r = 20, b = 15, l = 20, unit = "pt")
  )

print(p)
out <- "../hbond-v2/H-bond-E7wt1-E7apl1-v2"

# PDF：矢量，cairo_pdf 能正确嵌入字体和 ± / Å 等字符
ggsave(paste0(out, ".pdf"), plot = p, width = 8.5, height = 6, device = cairo_pdf)

# PNG：600 dpi
ggsave(paste0(out, ".png"), plot = p, width = 8.5, height = 6, dpi = 600)

# TIFF：600 dpi + LZW 压缩（投稿主用）
ggsave(paste0(out, ".tif"), plot = p, width = 8.5, height = 6, dpi = 600,
       device = "tiff", compression = "lzw")