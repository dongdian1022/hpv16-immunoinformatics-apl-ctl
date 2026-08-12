theme_publication <- function() {
  ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14, colour = "black"),
      axis.text = ggplot2::element_text(size = 12, colour = "black"),
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
      strip.text = ggplot2::element_text(size = 14, face = "bold"),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 11),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.6),
      plot.margin = ggplot2::margin(10, 12, 10, 12)
    )
}