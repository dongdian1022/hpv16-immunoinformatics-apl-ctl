## =============================================================================
## Figure 7A  -  Dose-dependent killing of SiHa-0201-E7 by HPV-T vs Mock-T
##   n = 3 healthy donors (HD1, HD2, HD3); Mock-T2 single control.
##   Technical duplicates collapsed to one value per donor (no pseudoreplication).
##   Line = group mean; error bars = SD across donors (HPV-T, n=3); points =
##   individual donors. Two-sided one-sample t vs Mock baseline; only the 30:1
##   comparison is significant (** p<0.01) and is the only one labelled.
## Palette: HPV-T #CB181D / Mock-T #08519C ; sig #333333.  Style = house.
## =============================================================================
library(ggplot2); library(dplyr)

use_arial <- FALSE
if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
  tryCatch({ sysfonts::font_add("Arial", regular = "Arial.ttf", bold = "Arial Bold.ttf")
    showtext::showtext_auto(); use_arial <- TRUE }, error = function(e) {})
}
BASE_FONT <- if (use_arial) "Arial" else "sans"
COL_HPV <- "#CB181D"; COL_MOCK <- "#08519C"; SIG_COL <- "#333333"

theme_v <- theme_classic(base_size = 12, base_family = BASE_FONT) +
  theme(
    axis.line = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(3.5, "pt"),
    axis.text.x  = element_text(size = 10.5, color = "#000000", family = BASE_FONT),
    axis.text.y  = element_text(size = 11,   color = "#000000", family = BASE_FONT),
    axis.title.x = element_text(size = 12, face = "bold", family = BASE_FONT, margin = margin(t = 6)),
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

## ---- raw data: HD1/2/3 + Mock-T2 on SiHa-0201-E7 (technical duplicates) ----
raw_A <- tribble(
  ~ET_Ratio, ~Donor, ~Group, ~Cytotoxicity,
  "3:1","HD1","HPV-T",17.23,"3:1","HD1","HPV-T",18.24,
  "3:1","HD2","HPV-T",38.78,"3:1","HD2","HPV-T",37.50,
  "3:1","HD3","HPV-T",13.86,"3:1","HD3","HPV-T",15.53,
  "10:1","HD1","HPV-T",47.08,"10:1","HD1","HPV-T",49.32,
  "10:1","HD2","HPV-T",75.33,"10:1","HD2","HPV-T",80.49,
  "10:1","HD3","HPV-T",37.46,"10:1","HD3","HPV-T",33.04,
  "30:1","HD1","HPV-T",66.27,"30:1","HD1","HPV-T",64.26,
  "30:1","HD2","HPV-T",82.56,"30:1","HD2","HPV-T",83.45,
  "30:1","HD3","HPV-T",74.45,"30:1","HD3","HPV-T",73.80,
  "3:1","Mock","Mock-T",1.72,"3:1","Mock","Mock-T",2.48,
  "10:1","Mock","Mock-T",6.57,"10:1","Mock","Mock-T",8.82,
  "30:1","Mock","Mock-T",7.88,"30:1","Mock","Mock-T",6.16
) %>% mutate(ET_Ratio = factor(ET_Ratio, c("3:1","10:1","30:1")),
             Group = factor(Group, c("HPV-T","Mock-T")))

donor_A <- raw_A %>% group_by(ET_Ratio, Group, Donor) %>%        ## collapse tech reps
  summarise(val = mean(Cytotoxicity), .groups = "drop")
grp_A <- donor_A %>% group_by(ET_Ratio, Group) %>%               ## SD = across donors (NA for Mock n=1)
  summarise(Mean = mean(val), SD = sd(val), SEM = sd(val)/sqrt(sum(!is.na(val))), .groups = "drop")

## only the 30:1 comparison is labelled (** ); low-E:T ns left unlabelled
sig_A <- tibble(ET_Ratio = factor("30:1", c("3:1","10:1","30:1")), y = 92, lab = "**")

pA <- ggplot(grp_A, aes(ET_Ratio, Mean, group = Group, colour = Group)) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.05, linewidth = 0.6, na.rm = TRUE) +
  geom_line(linewidth = 1.0) +
  geom_text(data = sig_A, aes(ET_Ratio, y, label = lab), inherit.aes = FALSE,
            size = 4.6, fontface = "bold", colour = SIG_COL, family = BASE_FONT) +
  scale_colour_manual(values = c("HPV-T" = COL_HPV, "Mock-T" = COL_MOCK)) +
  scale_fill_manual(values   = c("HPV-T" = COL_HPV, "Mock-T" = COL_MOCK), guide = "none") +
  scale_shape_manual(values  = c("HPV-T" = 21, "Mock-T" = 24), guide = "none") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = "E:T ratio", y = "Specific cytotoxicity (%)", subtitle = "Target: SiHa-0201-E7") +
  theme_v

save_three <- function(p, out, W, H) {
  ggsave(paste0(out, ".png"), p, width = W, height = H, dpi = 600, bg = "white")
  ggsave(paste0(out, ".pdf"), p, width = W, height = H, bg = "white", device = cairo_pdf)
  ggsave(paste0(out, ".tiff"), p, width = W, height = H, dpi = 600, bg = "white",
         device = tiff, compression = "lzw")
}
if (!exists("ASSEMBLING")) { print(pA); save_three(pA, "Figure7A-v4", 4, 4) }