## =============================================================================
## Figure 7C  -  HLA class I dependence (W6/32 blockade). E:T = 30:1, 21 h.
##   Groups: Mock-T2 (single control) and HPV-T (pooled n = 3 donors HD1/2/3).
##   Technical duplicates collapsed per donor. Bars = mean; error bars = SD
##   (HPV-T across donors, n=3); points = individual donors (HPV-T) / wells (Mock).
##   HPV-T without vs W6/32 blocked: PAIRED t-test, ** (p = 0.006); ~78-94%
##   reduction per donor. Mock-T2 (n=1) shown for reference, not tested.
##   alpha = condition; single condition legend. Style = house. Cap ***.
## Palette: HPV-T #CB181D / Mock-T #08519C ; sig/points #333333
## =============================================================================
library(ggplot2); library(dplyr)

use_arial <- FALSE
if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
  tryCatch({ sysfonts::font_add("Arial", regular = "Arial.ttf", bold = "Arial Bold.ttf")
    showtext::showtext_auto(); use_arial <- TRUE }, error = function(e) {})
}
BASE_FONT <- if (use_arial) "Arial" else "sans"
COL_HPV <- "#CB181D"; COL_MOCK <- "#08519C"; PT_COL <- "#333333"

theme_v <- theme_classic(base_size = 12, base_family = BASE_FONT) +
  theme(
    axis.line = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(3.5, "pt"),
    axis.text.x  = element_text(size = 10.5, color = "#000000", family = BASE_FONT),
    axis.text.y  = element_text(size = 11,   color = "#000000", family = BASE_FONT),
    axis.title.y = element_text(size = 12, face = "bold", family = BASE_FONT, margin = margin(r = 8)),
    plot.subtitle = element_text(size = 10.5, face = "bold", hjust = 0.5, family = BASE_FONT, margin = margin(b = 3)),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 18, 8, 8),
    legend.position = c(0.02, 1.00), legend.justification = c(0, 1),
    legend.direction = "vertical", legend.title = element_blank(),
    legend.text = element_text(size = 9.5, family = BASE_FONT),
    legend.background = element_blank(), legend.key.size = unit(11, "pt"),
    legend.key.width = unit(18, "pt"), legend.spacing.y = unit(0.02, "cm"))

## ---- donor-level values (HPV-T: HD1/2/3 means; Mock-T2: two technical wells) ----
raw_C <- tribble(
  ~Group, ~Condition, ~val,
  "HPV-T","Without block",65.265,"HPV-T","Without block",83.005,"HPV-T","Without block",74.125,
  "HPV-T","W6/32 blocked",12.215,"HPV-T","W6/32 blocked",18.165,"HPV-T","W6/32 blocked",4.810,
  "Mock-T","Without block",7.88,"Mock-T","Without block",6.16,
  "Mock-T","W6/32 blocked",3.32,"Mock-T","W6/32 blocked",2.36
) %>% mutate(Group = factor(Group, c("Mock-T","HPV-T")),
             Condition = factor(Condition, c("Without block","W6/32 blocked")))

summary_C <- raw_C %>% group_by(Group, Condition) %>%
  summarise(Mean = mean(val), SD = sd(val), SEM = sd(val)/sqrt(sum(!is.na(val))), .groups = "drop")

dge <- position_dodge(0.8)
## significance bracket over the HPV-T pair only (group x = 2; bars at 1.8 / 2.2)
brk <- data.frame(x = c(1.8, 1.8, 2.2, 2.2), y = c(88, 91, 91, 88))

pC <- ggplot(summary_C, aes(Group, Mean, fill = Group, alpha = Condition)) +
  geom_col(position = dge, width = 0.6, color = NA) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM, group = interaction(Group, Condition)),
                position = dge, width = 0.10, linewidth = 0.6, color = PT_COL,
                alpha = 1, na.rm = TRUE, show.legend = FALSE) +
  geom_path(data = brk, aes(x, y), inherit.aes = FALSE, colour = PT_COL, linewidth = 0.5) +
  annotate("text", x = 2, y = 94, label = "**", size = 4.6, fontface = "bold",
           colour = PT_COL, family = BASE_FONT) +
  scale_fill_manual(values = c("Mock-T" = COL_MOCK, "HPV-T" = COL_HPV), guide = "none") +
  scale_alpha_manual(values = c("Without block" = 1.0, "W6/32 blocked" = 0.4),
                     guide = guide_legend(override.aes = list(fill = "grey20"))) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = NULL, y = "Specific cytotoxicity (%)", subtitle = "Target: SiHa-0201-E7  (E:T = 30:1)") +
  theme_v

save_three <- function(p, out, W, H) {
  ggsave(paste0(out, ".png"), p, width = W, height = H, dpi = 600, bg = "white")
  ggsave(paste0(out, ".pdf"), p, width = W, height = H, bg = "white", device = cairo_pdf)
  ggsave(paste0(out, ".tiff"), p, width = W, height = H, dpi = 600, bg = "white",
         device = tiff, compression = "lzw")
}
if (!exists("ASSEMBLING")) { print(pC); save_three(pC, "Figure7C-v4", 4, 4) }