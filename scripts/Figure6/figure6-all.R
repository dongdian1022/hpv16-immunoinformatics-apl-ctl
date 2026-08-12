## =============================================================================
## assemble_Figure6.R  -  Composite Figure 6 (1 x 3), AXIS-ALIGNED via patchwork.
##   Why patchwork (not image stitching): patchwork aligns the plot PANELS across
##   sub-plots, so the y-axes / 0 / 100 baselines and panel sizes line up even
##   though the bar panel (C) has no x-axis title while the line panels (A,B) do.
##   Run this after the three panel scripts are in the same folder.
## =============================================================================
library(ggplot2)
library(patchwork)

ASSEMBLING <- TRUE                 # suppresses each panel's standalone export
here <- "."                        # set to the folder holding the 3 scripts
source(file.path(here, "Figure7A-v4.R"))            # -> pA
source(file.path(here, "Figure7B-v4.R")) # -> pB
source(file.path(here, "Figure7C-v4.R"))           # -> pC  (BASE_FONT now defined)

figure6 <- (pA | pB | pC) +
  plot_layout(widths = c(1, 1, 1)) +        # equal panel widths
  plot_annotation(tag_levels = "A") &       # bold A / B / C
  theme(plot.tag = element_text(size = 15, face = "bold", family = BASE_FONT),
        plot.tag.position = c(0.02, 0.98))

W <- 12; H <- 4
ggsave("Figure6.pdf",  figure6, width = W, height = H, bg = "white", device = cairo_pdf)
ggsave("Figure6.png",  figure6, width = W, height = H, dpi = 600, bg = "white")
ggsave("Figure6.tiff", figure6, width = W, height = H, dpi = 600, bg = "white",
       device = tiff, compression = "lzw")

message("Figure6-v4.{pdf,png,tiff} written — panels axis-aligned by patchwork.")