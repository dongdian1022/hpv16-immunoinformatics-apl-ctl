## =============================================================================
## Figure 7B  -  Target specificity: HPV-T vs Mock-T on parental SiHa
##   n = 3 healthy donors (HD1, HD2, HD3); Mock-T single control.
##   DESCRIPTIVE ONLY (no asterisks): on parental SiHa both HPV-T and Mock-T
##   show only background lysis (<10%). Target specificity is quantified in the
##   TEXT by a paired cross-target test (same donors, SiHa-0201-E7 vs parental):
##   significant at 30:1 (paired t-test, p = 0.003); 10:1 / 3:1 ns.
##   Y kept 0-100 to share scale with 7A. Style = house.
## Palette: HPV-T #CB181D / Mock-T #08519C
## =============================================================================
library(ggplot2); library(dplyr)

use_arial <- FALSE
if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
  tryCatch({ sysfonts::font_add("Arial", regular = "Arial.ttf", bold = "Arial Bold.ttf")
    showtext::showtext_auto(); use_arial <- TRUE }, error = function(e) {})
}
BASE_FONT <- if (use_arial) "Arial" else "sans"
COL_HPV <- "#CB181D"; COL_MOCK <- "#08519C"
YMAX <- 100

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

## ---- raw data: HD1/2/3 + Mock-T on parental SiHa (technical duplicates) ----
raw_B <- tribble(
  ~ET_Ratio, ~Donor, ~Group, ~Cytotoxicity,
  "3:1","HD1","HPV-T",1.70,"3:1","HD1","HPV-T",1.51,
  "3:1","HD2","HPV-T",0.60,"3:1","HD2","HPV-T",1.01,
  "3:1","HD3","HPV-T",0.96,"3:1","HD3","HPV-T",1.39,
  "10:1","HD1","HPV-T",2.25,"10:1","HD1","HPV-T",2.30,
  "10:1","HD2","HPV-T",2.51,"10:1","HD2","HPV-T",2.87,
  "10:1","HD3","HPV-T",1.95,"10:1","HD3","HPV-T",1.77,
  "30:1","HD1","HPV-T",1.66,"30:1","HD1","HPV-T",3.31,
  "30:1","HD2","HPV-T",6.92,"30:1","HD2","HPV-T",8.69,
  "30:1","HD3","HPV-T",3.14,"30:1","HD3","HPV-T",4.63,
  "3:1","Mock","Mock-T",0.54,"3:1","Mock","Mock-T",9.18,
  "10:1","Mock","Mock-T",3.60,"10:1","Mock","Mock-T",3.98,
  "30:1","Mock","Mock-T",7.94,"30:1","Mock","Mock-T",7.82
) %>% mutate(ET_Ratio = factor(ET_Ratio, c("3:1","10:1","30:1")),
             Group = factor(Group, c("HPV-T","Mock-T")))

donor_B <- raw_B %>% group_by(ET_Ratio, Group, Donor) %>%
  summarise(val = mean(Cytotoxicity), .groups = "drop")
grp_B <- donor_B %>% group_by(ET_Ratio, Group) %>%
  summarise(Mean = mean(val), SD = sd(val), SEM = sd(val)/sqrt(sum(!is.na(val))), .groups = "drop")

pB <- ggplot(grp_B, aes(ET_Ratio, Mean, group = Group, colour = Group)) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.05, linewidth = 0.6, na.rm = TRUE) +
  geom_line(linewidth = 1.0) +
  scale_colour_manual(values = c("HPV-T" = COL_HPV, "Mock-T" = COL_MOCK)) +
  scale_fill_manual(values   = c("HPV-T" = COL_HPV, "Mock-T" = COL_MOCK), guide = "none") +
  scale_shape_manual(values  = c("HPV-T" = 21, "Mock-T" = 24), guide = "none") +
  scale_y_continuous(limits = c(0, YMAX), breaks = seq(0, YMAX, 20), expand = c(0, 0)) +
  labs(x = "E:T ratio", y = "Specific cytotoxicity (%)", subtitle = "Target: Parental SiHa") +
  theme_v

save_three <- function(p, out, W, H) {
  ggsave(paste0(out, ".png"), p, width = W, height = H, dpi = 600, bg = "white")
  ggsave(paste0(out, ".pdf"), p, width = W, height = H, bg = "white", device = cairo_pdf)
  ggsave(paste0(out, ".tiff"), p, width = W, height = H, dpi = 600, bg = "white",
         device = tiff, compression = "lzw")
}
if (!exists("ASSEMBLING")) { print(pB); save_three(pB, "Figure7B-v4", 4, 4) }